{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds, KindSignatures, ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
module Data.Matrix.Sparse.Static where

import Control.Monad.ST
import Data.Proxy
import Data.Vector.Unboxed (Vector)
import qualified Data.Matrix.Sparse           as M
import qualified Numeric.LinearAlgebra.Static as Dense
import GHC.TypeLits

import qualified Data.Vector.Unboxed         as V
import qualified Data.Vector.Unboxed.Mutable as MV
import Data.STRef.Strict
import Data.MonoTraversable (Element)
import qualified Data.Vector.Sparse.Static as SV
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as MVS
import qualified Data.Vector.Generic as VG
import Control.Monad (forM_)

newtype Matrix (n :: Nat) (p :: Nat) a = Matrix (M.Matrix V.Vector a)
  deriving (Eq, Show, Read)

type instance Element (Matrix n p a) = a

{-# INLINE sparsify #-}
sparsify
  :: forall (n :: Nat) (p :: Nat).
     (KnownNat n, KnownNat p)
  => ((Int, Int) -> Double -> Maybe Double)
  -> Dense.L n p
  -> Matrix n p Double
sparsify f mat = Matrix (runST go)

  where maxLen = n * p
        n = fromIntegral $ natVal (Proxy :: Proxy n)
        p = fromIntegral $ natVal (Proxy :: Proxy p)
        go :: forall s. ST s (M.Matrix Vector Double)
        go = do
          rowIds <- MV.unsafeNew maxLen
          colIds <- MV.unsafeNew maxLen
          vals   <- MV.unsafeNew maxLen
          cntRef <- newSTRef 0
          flip Dense.imapLM_ mat $ \(i, j) v ->
            case f (i, j) v of
              Nothing -> return ()
              Just v' -> do
                cnt <- readSTRef cntRef
                writeSTRef cntRef (cnt+1)
                MV.unsafeWrite rowIds cnt i
                MV.unsafeWrite colIds cnt j
                MV.unsafeWrite vals   cnt v'
          cnt <- readSTRef cntRef
          rowIds' <- V.unsafeFreeze (MV.unsafeSlice 0 cnt rowIds)
          colIds' <- V.unsafeFreeze (MV.unsafeSlice 0 cnt colIds)
          vals'   <- V.unsafeFreeze (MV.unsafeSlice 0 cnt vals)
          return $ M.compress n p rowIds' colIds' vals'

{-# INLINE fromList #-}
fromList
  :: forall (n ::Nat) (p :: Nat) a.
     (KnownNat n, KnownNat p, Num a, MV.Unbox a)
  => [(Int, Int, a)] -> Matrix n p a
fromList xs = Matrix $ M.compress n p rowIds colIds vals
  where (rs, cs, vs) = unzip3 xs
        rowIds = V.fromList rs
        colIds = V.fromList cs
        vals   = V.fromList vs
        n = fromIntegral $ natVal (Proxy :: Proxy n)
        p = fromIntegral $ natVal (Proxy :: Proxy p)

{-# INLINE toList #-}
toList
  :: (KnownNat n, KnownNat p, MV.Unbox a)
  => Matrix n p a -> [(Int, Int, a)]
toList (Matrix m) =
  go (V.toList $ M.pointers m)
     (V.toList $ M.indices m)
     (V.toList $ M.values m)
     0

  where go _   [] _ _ = []
        go []  _  _ _ = []
        go [_] _  _ _ = []
        go (x:y:zs) rs vals j
          | x == y = go (y:zs) rs vals (j+1)
          | otherwise =
              let colLen = y - x
                  (colRowIds, otherRowIds) = splitAt colLen rs
                  (colVals, otherVals)     = splitAt colLen vals
                  xs = zipWith (\r v -> (r, j, v)) colRowIds colVals
              in xs ++ go (y:zs) otherRowIds otherVals (j+1)

{-# INLINE imap #-}
imap
  :: forall (n :: Nat) (p :: Nat) a b.
     (MV.Unbox a, MV.Unbox b)
  => (Int -> Int -> a -> b)
  -> Matrix n p a
  -> Matrix n p b
imap f (Matrix m) = Matrix $ m { M.values = values' }

  where values' = runST (go 0 0 =<< MV.unsafeNew n :: ST s (Vector b))
        n = V.length vs
        vs = M.values m
        rowIds = M.indices m
        colIds = M.pointers m
        go :: forall s. Int -> Int -> MV.MVector s b -> ST s (Vector b)
        go !col !valIdx mv
          -- we've processed all values (and/or all columns)
           | valIdx >= n = V.unsafeFreeze mv
          -- no nnz element in this column
           | (colIds V.! col) == (colIds V.! (col+1)) = go (col+1) valIdx mv
          -- otherwise, we've got some work to do
           | otherwise = do
              let !elemsInCol = (colIds V.! (col+1)) - (colIds V.! col)
              forM_ [0..(elemsInCol-1)] $ \i -> do
                let row = rowIds V.! (valIdx+i)
                    v   = vs V.! (valIdx+i)
                MV.unsafeWrite mv (valIdx+i) (f row col v)
              go (col+1) (valIdx+elemsInCol) mv

{-# INLINE add #-}
add :: (Num a, MV.Unbox a) => Matrix n p a -> Matrix n p a -> Matrix n p a
add (Matrix a) (Matrix b) = Matrix (a+b)

{-# INLINE mul #-}
mul
  :: (Num a, MV.Unbox a)
  => Matrix n p a -> Matrix p q a -> Matrix n q a
mul (Matrix a) (Matrix b) = Matrix (a * b)

nonZeros :: MV.Unbox a => Matrix n p a -> Int
nonZeros (Matrix m) = V.length (M.values m)

ident :: forall (n :: Nat) a. (KnownNat n, Num a, MV.Unbox a) => Matrix n n a
ident = Matrix $ M.Matrix
  { M.nrows = n
  , M.ncols = n
  , M.pointers = V.fromList $ [0..n] -- [0..(n-1)] ++ [n] to indicate number of nnz
  , M.indices  = V.fromList [0..(n-1)]
  , M.values   = V.replicate n 1
  }

  where n = fromIntegral $ natVal (Proxy :: Proxy n)

extractCol :: MV.Unbox a => Matrix n p a -> Int -> SV.V n a
extractCol (Matrix m) j = SV.V (M.unsafeSlice m j)

extractDenseCol :: Matrix n p Double -> Int -> Dense.R n
extractDenseCol m j = SV.toDense (extractCol m j)

nnzCol :: Matrix n p a -> Int -> Int
nnzCol (Matrix m) j = (ptrs V.! (j+1)) - (ptrs V.! j)
  where ptrs = M.pointers m

{-# INLINE transpose #-}
transpose :: MV.Unbox a => Matrix n p a -> Matrix p n a
transpose (Matrix m) = Matrix (M.transpose m)

{-# INLINE getCols #-}
getCols
  :: forall (n :: Nat) (p :: Nat) a.
     (KnownNat p, MV.Unbox a)
  => Matrix n p a -> [SV.V n a]
getCols m = map (extractCol m) [0..(p-1)]
  where p = fromIntegral $ natVal (Proxy :: Proxy p)

{-# INLINE getRows #-}
getRows
  :: forall (n :: Nat) (p :: Nat) a.
     (KnownNat n, MV.Unbox a)
  => Matrix n p a -> [SV.V p a]
getRows m = getCols (transpose m)

{-# INLINE sum #-}
sum :: (Num a, MV.Unbox a) => Matrix n p a -> a
sum (Matrix m) = V.sum (M.values m)

{-# INLINE zeroing #-}
zeroing
  :: (Num a, MV.Unbox a)
  => Int -> Int -> Matrix n p a -> (Matrix n p a -> r) -> r
zeroing i j mat@(Matrix m) f
  | Just k <- V.findIndex (==i) s_ids = runST $ do
      mvals <- V.unsafeThaw vals
      val <- MV.unsafeRead mvals (beg+k)
      MV.unsafeWrite mvals (beg+k) 0
      vals' <- V.unsafeFreeze mvals
      case f (Matrix $ m { M.values = vals' }) of
        r -> MV.unsafeWrite mvals (beg+k) val >> return r
  | otherwise = f mat

  where ptrs = M.pointers m
        ids  = M.indices m
        vals = M.values m
        beg  = ptrs V.! j
        end  = ptrs V.! (j+1)
        s_ids = V.unsafeSlice beg (end - beg) ids

{-# INLINE asRow #-}
asRow
  :: forall (n :: Nat) a. (KnownNat n, Num a, MV.Unbox a)
  => SV.V n a -> Matrix 1 n a
asRow v = Matrix $ M.compress 1 n rows cols vals

  where n  = SV.dim v
        vals = SV.values v
        nnz = SV.nnz v
        rows = V.replicate nnz 0
        cols = SV.indices v

asColumn :: forall (n :: Nat) a. KnownNat n => SV.V n a -> Matrix n 1 a
asColumn v = Matrix $ M.Matrix 1 n ptrs ids vs
  where ids = SV.indices v
        nnz = V.length ids
        vs  = SV.values v
        ptrs = V.fromList [0, nnz]
        n   = fromIntegral $ natVal (Proxy :: Proxy n)

dia :: forall (n :: Nat). KnownNat n => Dense.R n -> Matrix n n Double
dia (Dense.R (Dense.Dim v)) = Matrix $ M.Matrix n n ptrs ids vs
  where vs = VG.convert v
        ptrs = V.enumFromTo 0 n
        ids = V.enumFromTo 0 n
        n = fromIntegral $ natVal (Proxy :: Proxy n)

{-# INLINE at #-}
at :: (MV.Unbox a, Num a) => Matrix n p a -> (Int, Int) -> a
at m (i, j) = extractCol m j SV.! i

{-# INLINE densify #-}
densify
  :: forall (n :: Nat) (p :: Nat).
     (KnownNat n, KnownNat p)
  => Matrix n p Double -> Dense.L n p
densify (Matrix m) = Dense.L . Dense.Dim . Dense.Dim $ M.pack m

-- symmetric

-- we just store the (strict) upper triangle (col > row) sparsely,
-- + diagonal densely

data Sym n = Sym
  { diag :: !(Dense.R n)
  , uptr :: !(Matrix n n Double)
  } deriving Show

{-# INLINE fromListSym #-}
-- the 2nd arg should only include (i, j) pairs such that i < j
fromListSym :: KnownNat n => [Double] -> [(Int, Int, Double)] -> Sym n
fromListSym d rest = Sym d' u
  where d' = Dense.fromList d
        u = fromList rest

{-# INLINE imapSym #-}
imapSym :: (Int -> Int -> Double -> Double) -> Sym n -> Sym n
imapSym f (Sym d u) = Sym (mapDiag d) (imap f u)
  where mapDiag (Dense.R (Dense.Dim v)) = Dense.R $ Dense.Dim $
          VS.imap (\i a -> f i i a) v

{-# INLINE symNnzCol #-}
symNnzCol :: Sym n -> Int -> Int
symNnzCol (Sym (Dense.R (Dense.Dim d)) u) j =
  (if d VS.! j /= 0 then 1 else 0) +
  (nnzCol u j + nnzRowAfter u j j)

-- how many non zero entries do we have in 'm' at row 'row',
-- starting the count _after_ column 'col'
{-# INLINE nnzRowAfter #-}
nnzRowAfter :: Matrix n n a -> Int -> Int -> Int
nnzRowAfter (Matrix m) row col = V.foldl' cnt 0 v

  where startIdx = M.pointers m V.! (col+1)
        len      = V.length (M.indices m)
        v        = V.unsafeSlice startIdx (len - startIdx) (M.indices m)
        cnt k rownum = if rownum == row then k+1 else k

{-# INLINE symAt #-}
symAt :: Sym n -> (Int, Int) -> Double
symAt m (i, j)
  | i < j = at (uptr m) (i, j)
  | i == j = case diag m of
      Dense.R (Dense.Dim v) -> v VS.! i
  | otherwise = symAt m (j, i)

{-# INLINE symNnz #-}
symNnz :: Sym n -> Int
symNnz m = nonZeros (uptr m) * 2 + VS.foldl' cnt 0 (unR $ diag m)
  where cnt k v = if v /= 0 then k+1 else k
        unR (Dense.R (Dense.Dim v)) = v

{-# INLINE zeroingV #-}
zeroingV :: Int -> Dense.R n -> (Dense.R n -> r) -> r
zeroingV i (Dense.R (Dense.Dim v)) f = runST $ do
  mv <- VS.unsafeThaw v
  a  <- MVS.unsafeRead mv i
  MVS.unsafeWrite mv i 0
  v' <- VS.unsafeFreeze mv
  case f (Dense.R (Dense.Dim v')) of
    r -> MVS.unsafeWrite mv i a >> return r

{-# INLINE zeroingSym #-}
zeroingSym :: Int -> Int -> Sym n -> (Sym n -> r) -> r
zeroingSym i j m@(Sym d u) f
  | i < j = zeroing i j u $ \uptr' -> f (Sym d uptr')
  | i == j = zeroingV i d $ \d' -> f (Sym d' u)
  | otherwise = zeroingSym j i m f

-- implements SssSpMVSerial from:
--   https://kkourt.io/papers/ipdps13.pdf
{-# INLINE mulVSym #-}
mulVSym :: forall (n :: Nat). KnownNat n => Sym n -> Dense.R n -> Dense.R n
mulVSym (Sym (Dense.R (Dense.Dim d)) (Matrix u)) (Dense.R (Dense.Dim v)) = runST $ do
  out <- MVS.new n
  forM_ [0..(n-1)] $ \r -> do
    MVS.unsafeWrite out r $ (d VS.! r) * (v VS.! r)
    forM_ [ (ptrs V.! r) .. ((ptrs V.! (r+1)) - 1) ] $ \j -> do
      let c = ids V.! j
      MVS.unsafeModify out (\a -> a + (vs V.! j) * (v VS.! c)) r
      MVS.unsafeModify out (\a -> a + (vs V.! j) * (v VS.! r)) c
  Dense.R . Dense.Dim <$> VS.unsafeFreeze out

  where n = fromIntegral $ natVal (Proxy :: Proxy n)
        ptrs = M.pointers u
        ids  = M.indices u
        vs   = M.values u

{-# INLINE unSym #-}
unSym :: KnownNat n => Sym n -> Matrix n n Double
unSym (Sym d u) = dia d `add` u `add` transpose u

{-# INLINE toListSym #-}
toListSym
  :: KnownNat n => Sym n -> [(Int, Int, Double)]
toListSym (Sym (Dense.R (Dense.Dim d)) u) =
  [ (i, i, a) | (i, a) <- Prelude.zip [0..] (VS.toList d) ] ++
  concatMap (\(i, j, a) -> [(i, j, a), (j, i, a)]) (toList u)
