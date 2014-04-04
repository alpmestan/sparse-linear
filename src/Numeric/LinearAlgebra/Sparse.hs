{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE TypeFamilies #-}

module Numeric.LinearAlgebra.Sparse
    ( Matrix
    , OrderK(..), OrderR(..)
    , FormatK(..), FormatR(..)
    , slices, rows, cols, slice
    , pack, reorder
    , mulV, mulVM, mul, add
    ) where

import Control.Exception (assert)
import Control.Lens
import Control.Monad (liftM, liftM2, when)
import Control.Monad.Primitive (PrimMonad(..))
import Control.Monad.ST (runST)
import Data.AEq
import Data.List (foldl')
import Data.Maybe (fromMaybe)
import Data.Monoid ((<>))
import Data.Ord (comparing)
import qualified Data.Vector.Generic as V
import qualified Data.Vector.Generic.Mutable as MV
import Data.Vector.Unboxed (Unbox, Vector)
import qualified Data.Vector.Unboxed as U
import Data.Vector.Unboxed.Mutable (MVector)
import qualified Data.Vector.Unboxed.Mutable as MU
import Data.Vector.Algorithms.Intro (sortBy)
import Data.Vector.Algorithms.Search (binarySearchL)

import Data.Proxy.PolyKind

-- | Matrix order
data OrderK
    = Row -- ^ row-major
    | Col -- ^ column-major

-- | Matrix formats
data FormatK
    = U -- ^ uncompressed
    | C -- ^ compressed

class OrderR (ord :: OrderK) where
    -- | Extract the major index from a tuple in (row, column, ...) order.
    major :: (Field1 s s a a, Field2 s s a a)
          => Lens' (Tagged ord s) a
    -- | Extract the minor index from a tuple in (row, column, ...) order
    minor :: (Field1 s s a a, Field2 s s a a)
          => Lens' (Tagged ord s) a
    -- | Extract the row index from a tuple in (major, minor, ...) order.
    row :: (Field1 s s a a, Field2 s s a a)
        => Lens' (Tagged ord s) a
    -- | Extract the column index from a tuple in (major, minor, ...) order
    col :: (Field1 s s a a, Field2 s s a a)
        => Lens' (Tagged ord s) a

{-# INLINE comparator #-}
comparator  :: (Field1 s s a a, Field2 s s a a, Ord a, OrderR ord)
            => Proxy ord -> s -> s -> Ordering
comparator witness a' b' =
    let a = tag witness a'
        b = tag witness b'
    in comparing (view major) a b <> comparing (view minor) a b

instance OrderR Row where
    major f = fmap (tag Proxy) . _1 f . untag
    minor f = fmap (tag Proxy) . _2 f . untag
    row f = fmap (tag Proxy) . _1 f . untag
    col f = fmap (tag Proxy) . _2 f . untag
    {-# INLINE major #-}
    {-# INLINE minor #-}
    {-# INLINE row #-}
    {-# INLINE col #-}

instance OrderR Col where
    major f = fmap (tag Proxy) . _2 f . untag
    minor f = fmap (tag Proxy) . _1 f . untag
    row f = fmap (tag Proxy) . _2 f . untag
    col f = fmap (tag Proxy) . _1 f . untag
    {-# INLINE major #-}
    {-# INLINE minor #-}
    {-# INLINE row #-}
    {-# INLINE col #-}

-- | Compressed sparse format
data Cx a = Cx !Int -- ^ minor dimension
               !(Vector Int) -- ^ starting indices of each major slice
               !(Vector (Int, a)) -- ^ (minor index, coefficient)

instance (Eq a, Unbox a) => Eq (Cx a) where
    (==) (Cx mnrA ixsA valsA) (Cx mnrB ixsB valsB) =
        mnrA == mnrB && ixsA == ixsB && valsA == valsB

instance (AEq a, Unbox a) => AEq (Cx a) where
    (===) (Cx mnrA ixsA valsA) (Cx mnrB ixsB valsB) =
        let (nsA, coeffsA) = U.unzip valsA
            (nsB, coeffsB) = U.unzip valsB
        in mnrA == mnrB && ixsA == ixsB && nsA == nsB
            && U.and (U.zipWith (===) coeffsA coeffsB)

    (~==) (Cx mnrA ixsA valsA) (Cx mnrB ixsB valsB) =
        let (nsA, coeffsA) = U.unzip valsA
            (nsB, coeffsB) = U.unzip valsB
        in mnrA == mnrB && ixsA == ixsB && nsA == nsB
            && U.and (U.zipWith (~==) coeffsA coeffsB)

-- | Uncompressed sparse format
data Ux a = Ux !Int -- ^ row dimension
               !Int -- ^ column dimension
               !(Vector (Int, Int, a))
               -- ^ (row index, column index, coefficient)

instance (Eq a, Unbox a) => Eq (Ux a) where
    (==) (Ux rA cA valsA) (Ux rB cB valsB) =
        rA == rB && cA == cB && valsA == valsB

instance (AEq a, Unbox a) => AEq (Ux a) where
    (===) (Ux rA cA valsA) (Ux rB cB valsB) =
        let (rsA, csA, coeffsA) = U.unzip3 valsA
            (rsB, csB, coeffsB) = U.unzip3 valsB
        in rA == rB && cA == cB && rsA == rsB && csA == csB
            && U.and (U.zipWith (===) coeffsA coeffsB)

    (~==) (Ux rA cA valsA) (Ux rB cB valsB) =
        let (rsA, csA, coeffsA) = U.unzip3 valsA
            (rsB, csB, coeffsB) = U.unzip3 valsB
        in rA == rB && cA == cB && rsA == rsB && csA == csB
            && U.and (U.zipWith (~==) coeffsA coeffsB)

data family Matrix :: FormatK -> OrderK -> * -> *
newtype instance Matrix C ord a = MatC (Tagged ord (Cx a))
newtype instance Matrix U ord a = MatU (Tagged ord (Ux a))

class FormatR (fmt :: FormatK) where
    -- | The dimensions of a matrix in (row, column) order.
    dim :: OrderR ord => Matrix fmt ord a -> (Int, Int)

    -- | The dimensions of a matrix in format-specific (major, minor)
    -- order.
    dimF :: OrderR ord => Matrix fmt ord a -> (Int, Int)

    -- | The number of non-zero entries in the matrix.
    nonzero :: Unbox a => Matrix fmt ord a -> Int

    compress :: (OrderR ord, Unbox a) => Matrix fmt ord a -> Matrix C ord a

    decompress :: (OrderR ord, Unbox a) => Matrix fmt ord a -> Matrix U ord a

    transpose :: (OrderR ord, Unbox a) => Matrix fmt ord a -> Matrix fmt ord a

instance FormatR U where
    dim (MatU (Tagged (Ux r c _))) = (r, c)

    dimF mat@(MatU ux) =
      -- TODO: Find a lens-y way to do this
      let _dim = copyTag ux $ dim mat
          m = view major _dim
          n = view minor _dim
      in (m, n)

    nonzero (MatU ux) =
      let Ux _ _ vals = untag ux
      in U.length vals

    compress mat@(MatU ux) =
        MatC $ unproxy $ \witness ->
          let (Ux _ _ vals) = proxy ux witness
              (rows_, cols_, coeffs) = U.unzip3 vals
              (m, n) = dimF mat
              majors = view major $ tag witness (rows_, cols_)
              minors = view minor $ tag witness (rows_, cols_)
              ixs = runST $ do
                majorsM <- U.thaw majors
                U.generateM m $ binarySearchL majorsM
          in Cx n ixs $ U.zip minors coeffs

    decompress x = x

    transpose (MatU ux) =
        MatU $ unproxy $ \witness ->
            let (Ux r c vals) = proxy ux witness
            in sortUx witness $ Ux c r $ U.map (\(x, y, z) -> (y, x, z)) vals

    {-# INLINE dim #-}
    {-# INLINE dimF #-}
    {-# INLINE nonzero #-}
    {-# INLINE compress #-}
    {-# INLINE decompress #-}
    {-# INLINE transpose #-}

instance FormatR C where
    dimF (MatC (Tagged (Cx n ixs _))) = (U.length ixs, n)

    dim mat@(MatC cx) =
      -- TODO: Find a lens-y way to do this
      let _dim = copyTag cx $ dimF mat
          r = view row _dim
          c = view col _dim
      in (r, c)

    compress x = x

    nonzero (MatC cx) =
      let Cx _ _ vals = untag cx
      in U.length vals

    decompress mat@(MatC cx) =
        MatU $ unproxy $ \witness ->
          let (Cx _ majorIxs valsC) = proxy cx witness
              (minors, coeffs) = U.unzip valsC
              (m, _) = dimF mat
              majorLengths =
                U.generate m $ \r ->
                  let this = majorIxs U.! r
                      next = fromMaybe (U.length valsC) (majorIxs U.!? (succ r))
                  in next - this
              majors =
                U.concatMap (\(r, len) -> U.replicate len r)
                $ U.indexed majorLengths
              rows_ = view row $ tag witness (majors, minors)
              cols_ = view col $ tag witness (majors, minors)
              (nr, nc) = dim mat
          in Ux nr nc $ U.zip3 rows_ cols_ coeffs

    transpose = compress . transpose . decompress

    {-# INLINE dim #-}
    {-# INLINE dimF #-}
    {-# INLINE nonzero #-}
    {-# INLINE compress #-}
    {-# INLINE decompress #-}
    {-# INLINE transpose #-}

instance (Eq a, Unbox a) => Eq (Matrix C ord a) where
    (==) (MatC a) (MatC b) = untag a == untag b

instance (Eq a, Unbox a) => Eq (Matrix U ord a) where
    (==) (MatU a) (MatU b) = untag a == untag b

instance (AEq a, Unbox a) => AEq (Matrix C ord a) where
    (===) (MatC a) (MatC b) = untag a === untag b
    (~==) (MatC a) (MatC b) = untag a ~== untag b

instance (AEq a, Unbox a) => AEq (Matrix U ord a) where
    (===) (MatU a) (MatU b) = untag a === untag b
    (~==) (MatU a) (MatU b) = untag a ~== untag b

{-# INLINE generate #-}
generate :: Int -> (Int -> a) -> [a]
generate len f = map f $ take len $ [0..]

-- | I was going to put this in a class to abstract over the matrix format,
-- but I realized that I would just end up compressing the matrix before
-- slicing, anyway. It's a Fold, not a Traversable, because it's not
-- suitable for modifying the matrix. If you want to do that, you'll
-- probably get better speed by manipulating the structure directly. This
-- is just for consumption.
{-# INLINE slices #-}
slices :: (Unbox a, OrderR ord) => Fold (Matrix C ord a) (Vector (Int, a))
slices = folding $ \mat -> generate (fst $ dimF mat) $ slice mat

{-# INLINE rows #-}
rows :: Unbox a => Fold (Matrix C Row a) (Vector (Int, a))
rows = slices

{-# INLINE cols #-}
cols :: Unbox a => Fold (Matrix C Col a) (Vector (Int, a))
cols = slices

{-# INLINE sortUx #-}
sortUx :: (OrderR ord, Unbox a) => Proxy ord -> Ux a -> Ux a
sortUx witness (Ux nr nc vals) =
    Ux nr nc $ U.modify (sortBy $ comparator witness) vals

{-# INLINE reorder #-}
reorder :: (OrderR ord', Unbox a) => Matrix U ord a -> Matrix U ord' a
reorder (MatU mat) =
    MatU $ unproxy $ \witness -> sortUx witness $ untag mat

{-# INLINE pack #-}
pack :: (OrderR ord, Unbox a)
     => Int -> Int -> Vector (Int, Int, a) -> Matrix U ord a
pack r c v = MatU $ unproxy $ \witness -> sortUx witness $ Ux r c v

{-# INLINE mulV #-}
mulV  :: (Num a, Unbox a, V.Vector v a) => Matrix C Row a -> v a -> v a
mulV mat xs_ =
    assert (c == U.length xs)
    $ V.convert $ U.create $ do
      ys <- MU.new r
      iforMOf_ (indexing rows) mat $ \ixR _row -> do
        let (cols_, coeffs) = U.unzip _row
        MU.write ys ixR
          $ U.sum $ U.zipWith (*) coeffs
          $ U.backpermute xs cols_
      return ys
  where
    xs = V.convert xs_
    (r, c) = dim mat

{-# INLINE mulVM #-}
mulVM :: (MV.MVector v a, Num a, PrimMonad m, Unbox a)
      => Matrix C Row a -> v (PrimState m) a -> v (PrimState m) a -> m ()
mulVM mat src dst =
    assert (c == MV.length src)
    $ assert (r == MV.length dst)
    $ iforMOf_ (indexing rows) mat $ \ixR _row -> do
      let (_cols, _coeffs) = U.unzip _row
      x <- liftM (U.sum . U.zipWith (*) _coeffs) $ U.mapM (MV.read src) _cols
      MV.write dst ixR x
  where
    (r, c) = dim mat

{-# INLINE mul #-}
mul :: (Num a, OrderR ord, Unbox a)
    => Matrix C Col a -> Matrix C Row a -> Matrix C ord a
mul a b =
    assert (inner == inner')
    $ foldl' add empty $ generate inner $ \i -> expand (slice a i) (slice b i)
  where
    (inner, left) = dimF a
    (inner', right) = dimF b
    empty = compress $ pack left right U.empty
    expand :: (Num a, OrderR ord, Unbox a)
           => Vector (Int, a) -> Vector (Int, a) -> Matrix C ord a
    expand ls rs = MatC $ unproxy $ \witness ->
      let dimC = tag witness (left, right)
          lsrs = tag witness (ls, rs)
          mjr = view major dimC
          mnr = view minor dimC
          ms = view major lsrs
          ns = view minor lsrs
          vals = U.concatMap (\(_, x) -> U.map (\(n, y) -> (n, x * y)) ns) ms
          stride = U.length ns
          ixs = U.iterateN mjr (+ stride) 0
      in Cx mnr ixs vals

{-# INLINE slice #-}
slice :: Unbox a => Matrix C ord a -> Int -> Vector (Int, a)
slice (MatC cx) i =
    let Cx _ starts vals = untag cx
        start = starts U.! i
        end = fromMaybe (U.length vals) $ starts U.!? i
    in assert (i < U.length starts)
    $ (U.slice start $ end - start) vals

{-# INLINE add #-}
add :: (Num a, OrderR ord, Unbox a)
    => Matrix C ord a -> Matrix C ord a -> Matrix C ord a
add a b =
    let (majorA, minorA) = dimF a
        (valsC, ixsC) = runST $ do
          vals <- MU.new $ nonzero a + nonzero b
          ixs <- MU.new majorA
          let go i start
                | i < majorA = do
                  MU.unsafeWrite ixs i start
                  let sliceA = slice a i
                      sliceB = slice b i
                  len <- addSlicesInto
                    (MU.slice start (U.length sliceA + U.length sliceB) vals)
                    sliceA sliceB
                  go (succ i) (start + len)
                | otherwise = return start
          len <- go 0 0
          vals_ <- U.unsafeFreeze $ MU.unsafeSlice 0 len vals
          ixs_ <- U.unsafeFreeze ixs
          return (vals_, ixs_)
    in assert (dimF a == dimF b)
        $ MatC $ tag Proxy $ Cx minorA ixsC valsC
  where
    addSlicesInto :: (Monad m, Num a, PrimMonad m, Unbox a)
                  => MVector (PrimState m) (Int, a)
                  -> Vector (Int, a)
                  -> Vector (Int, a)
                  -> m Int
    addSlicesInto dst src1 src2 = do
      let len1 = U.length src1
          len2 = U.length src2
          len = MU.length dst
      MU.move (MU.slice 0 len1 dst) =<< U.thaw src1
      MU.move (MU.slice len1 (len1 + len2) dst) =<< U.thaw src2
      sortBy (comparing fst) dst

      -- Accumulate elements in same minor dimension by adding their
      -- coefficients together. Set the minor dimension of duplicate
      -- elements to (-1) so we can find them later.
      let (ns, xs) = MU.unzip dst -- unzip to avoid looking up both fields
          accumulate update check =
            when (check < len) $ do
              -- do the elements have the same minor dimension?
              -- i.e., is one a duplicate?
              dup <- liftM2 (==) (MU.read ns update) (MU.read ns check)
              if dup
                then do
                  -- add the duplicate coefficients
                  x <- liftM2 (+) (MU.read xs update) (MU.read xs check)
                  -- insert new coefficient in place of first duplicate
                  MU.write xs update x
                  -- mark second duplicate for removal
                  MU.write ns check (-1)
                  accumulate update (succ check)
                else accumulate check (succ check)
      accumulate 0 1

      sortBy (comparing fst) dst
      start <- binarySearchL (fst $ MU.unzip dst) 0
      let len' = len - start
      when (start > 0) $ MU.move (MU.slice 0 len' dst) (MU.slice start len' dst)
      return len'