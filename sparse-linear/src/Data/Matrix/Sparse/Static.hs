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
import Control.Monad (forM_, when)
import qualified Data.Vector.Sparse.Static as SV
import Debug.Trace (trace, traceShow)

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
        vs_len = V.length vs
        rowIds = M.indices m
        colIds = M.pointers m
        colIds_len = V.length colIds
        go :: Int -> Int -> MV.MVector s b -> ST s (Vector b)
        go !col !valIdx mv
          -- we've processed all values (and/or all columns)
           | valIdx >= vs_len || col >= colIds_len = V.unsafeFreeze mv
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

{-# INLINE izipR #-}
-- an indexed-zip where the result will have the same structure as the 2nd
-- Matrix argument
izipR
  :: forall (n :: Nat) (p :: Nat) a b c.
     (MV.Unbox a, MV.Unbox b, MV.Unbox c)
  => a
  -> (Int -> Int -> a -> b -> c)
  -> Matrix n p a
  -> Matrix n p b
  -> Matrix n p c
izipR za f (Matrix m1) (Matrix m2) = Matrix $ M.Matrix
  { M.ncols = ncols
  , M.nrows = M.nrows m2
  , M.pointers = colIds2
  , M.indices = rowIds2
  , M.values = vs
  }

  where ncols = M.ncols m2
        nvals2 = V.length vals2
        rowIds2 = M.indices m2
        rowIds1 = M.indices m1
        colIds1 = M.pointers m1
        colIds2 = M.pointers m2
        vals2 = M.values m2
        vals1 = M.values m1
        vs = runST $ do
          mvs <- MV.unsafeNew nvals2
          go 0 0 0 mvs
          V.unsafeFreeze mvs

        go :: forall s.
              Int -> Int -> Int
           -> MV.MVector s c
           -> ST s ()
        go !col !valIdx1 !valIdx2 mv
          | (valIdx2 >= nvals2) = return ()
          | (colIds2 V.! col) == (colIds2 V.! (col+1)) = do
              let elemsInCol1 = (colIds1 V.! (col+1)) - (colIds1 V.! col)
              go (col+1) (valIdx1+elemsInCol1) valIdx2 mv
          | otherwise = do
              let !elemsInCol2 = (colIds2 V.! (col+1)) - (colIds2 V.! col)
                  !elemsInCol1 = (colIds1 V.! (col+1)) - (colIds1 V.! col)
                  goCol :: Int -> Int -> ST s ()
                  goCol !i1 !i2 = do
                    when (i2 < elemsInCol2) $ do
                        let !row2 = rowIds2 V.! (valIdx2 + i2)
                            !val2 = vals2 V.! (valIdx2 + i2)

                            midx1 = V.findIndex (>= row2)
                              (V.unsafeSlice (valIdx1+i1) elemsInCol1 rowIds1)

                            (!val1, !nextIdx1) = case midx1 of
                              Just i1'
                                | rowIds1 V.! (valIdx1+i1+i1') == row2 ->
                                    (vals1 V.! (valIdx1+i1+i1'), i1+i1'+1)
                                | otherwise ->
                                    (za, i1+i1')
                              Nothing -> (za, i1)

                        MV.unsafeWrite mv (valIdx2+i2) (f row2 col val1 val2)
                        goCol nextIdx1 (i2+1)
              goCol 0 0
              go (col+1) (valIdx1+elemsInCol1) (valIdx2+elemsInCol2) mv

{-# INLINE izipL #-}
-- an indexed-zip where the result will have the same structure as the 1st
-- Matrix argument
izipL
  :: forall (n :: Nat) (p :: Nat) a b c.
     (MV.Unbox a, MV.Unbox b, MV.Unbox c)
  => b
  -> (Int -> Int -> a -> b -> c)
  -> Matrix n p a
  -> Matrix n p b
  -> Matrix n p c
izipL zb f m1 m2 = izipR zb (\i j a b -> f i j b a) m2 m1

zipR
  :: (MV.Unbox a, MV.Unbox b, MV.Unbox c)
  => a
  -> (a -> b -> c)
  -> Matrix n p a
  -> Matrix n p b
  -> Matrix n p c
zipR za f m1 m2 = izipR za (\_ _ a b -> f a b) m1 m2

zipL
  :: (MV.Unbox a, MV.Unbox b, MV.Unbox c)
  => b
  -> (a -> b -> c)
  -> Matrix n p a
  -> Matrix n p b
  -> Matrix n p c
zipL zb f m1 m2 = izipL zb (\_ _ a b -> f a b) m1 m2

{-# INLINE izip #-}
-- an indexed-zip where the result will have the nonzero of structure
-- of the "union" of the two argument matrices (i.e an element will be
-- explicitly stored in the resulting matrix when either of the argument
-- matrices were explicitly storing the element at the given coordinate)

-- TODO: version that compresses on the fly?
--       with 'Int -> Int -> a -> b -> Maybe c' argument?
izip
  :: forall (n :: Nat) (p :: Nat) a b c.
     (MV.Unbox a, MV.Unbox b, MV.Unbox c)
  => a
  -> b
  -> (Int -> Int -> a -> b -> c)
  -> Matrix n p a
  -> Matrix n p b
  -> Matrix n p c
izip za zb f (Matrix m1) (Matrix m2) = Matrix $
   M.Matrix
  { M.ncols = ncols
  , M.nrows = nrows
  , M.pointers = cs
  , M.indices = rs
  , M.values = vs
  }

  where ncols = M.ncols m2
        nrows = M.nrows m2
        nvals1 = V.length vals1
        nvals2 = V.length vals2
        rowIds2 = M.indices m2
        rowIds1 = M.indices m1
        colIds1 = M.pointers m1
        colIds2 = M.pointers m2
        vals2 = M.values m2
        vals1 = M.values m1
        (rs, cs, vs) = runST $ do
          mcs <- MV.unsafeNew (ncols+1)
          mrs <- MV.unsafeNew (nvals1 + nvals2)
          mvs <- MV.unsafeNew (nvals1 + nvals2)
          MV.unsafeWrite mcs 0 0
          nnz <- go 0 0 0 0 mrs mcs mvs
          (,,) <$> V.unsafeFreeze (MV.unsafeSlice 0 nnz mrs)
               <*> V.unsafeFreeze mcs
               <*> V.unsafeFreeze (MV.unsafeSlice 0 nnz mvs)

        go :: forall s.
              Int -> Int -> Int -> Int
           -> MV.MVector s Int
           -> MV.MVector s Int
           -> MV.MVector s c
           -> ST s Int
        go !nnz !col !valIdx1 !valIdx2 mrs mcs mvs
          | (valIdx1 >= nvals1 && valIdx2 >= nvals2) = do
              forM_ [(col+1)..ncols] $ \i -> MV.unsafeWrite mcs i nnz
              return nnz

          | (colIds1 V.! col) == (colIds1 V.! (col+1)) &&
            (colIds2 V.! col) == (colIds2 V.! (col+1)) = do
              MV.unsafeWrite mcs (col+1) nnz
              go nnz (col+1) valIdx1 valIdx2 mrs mcs mvs

          | otherwise = do
              let !elemsInCol1 = (colIds1 V.! (col+1)) - (colIds1 V.! col)
                  !elemsInCol2 = (colIds2 V.! (col+1)) - (colIds2 V.! col)
                  rs1 = V.unsafeSlice valIdx1 elemsInCol1 rowIds1
                  vs1 = V.unsafeSlice valIdx1 elemsInCol1 vals1
                  rs2 = V.unsafeSlice valIdx2 elemsInCol2 rowIds2
                  vs2 = V.unsafeSlice valIdx2 elemsInCol2 vals2
                  goCol !i1 !i2 !nnzs
                    | i1 < elemsInCol1 && i2 < elemsInCol2 = do
                        let !r1 = rs1 V.! i1
                            !r2 = rs2 V.! i2
                        case compare r1 r2 of
                          LT -> do
                            MV.unsafeWrite mrs (nnz + nnzs) r1
                            MV.unsafeWrite mvs (nnz + nnzs) $ f r1 col (vs1 V.! i1) zb
                            goCol (i1+1) i2 (nnzs+1)
                          GT -> do
                            MV.unsafeWrite mrs (nnz + nnzs) r2
                            MV.unsafeWrite mvs (nnz + nnzs) $ f r2 col za (vs2 V.! i2)
                            goCol i1 (i2+1) (nnzs+1)
                          EQ -> do
                            MV.unsafeWrite mrs (nnz+nnzs) r1
                            MV.unsafeWrite mvs (nnz + nnzs) $ f r1 col (vs1 V.! i1) (vs2 V.! i2)
                            goCol (i1+1) (i2+1) (nnzs+1)
                    | i1 < elemsInCol1 = do
                        let remainingRows1 = V.unsafeSlice i1 (elemsInCol1 - i1) rs1
                            remainingVals1 = V.unsafeSlice i1 (elemsInCol1 - i1) vs1
                        V.izipWithM_ (\i r v -> MV.unsafeWrite mrs (nnz + nnzs + i) r
                                             >> MV.unsafeWrite mvs (nnz + nnzs + i) (f r col v zb)
                                     )
                                     remainingRows1 remainingVals1
                        return (nnzs + elemsInCol1 - i1)
                    | otherwise = do
                        let remainingRows2 = V.unsafeSlice i2 (elemsInCol2 - i2) rs2
                            remainingVals2 = V.unsafeSlice i2 (elemsInCol2 - i2) vs2
                        V.izipWithM_ (\i r v -> MV.unsafeWrite mrs (nnz + nnzs + i) r
                                             >> MV.unsafeWrite mvs (nnz + nnzs + i) (f r col za v)
                                     )
                                     remainingRows2 remainingVals2
                        return (nnzs + elemsInCol2 - i2)
              nnz' <- goCol 0 0 0
              MV.unsafeWrite mcs (col+1) (nnz + nnz')
              go (nnz + nnz') (col+1) (valIdx1+elemsInCol1) (valIdx2+elemsInCol2) mrs mcs mvs

zip
  :: (MV.Unbox a, MV.Unbox b, MV.Unbox c)
  => a
  -> b
  -> (a -> b -> c)
  -> Matrix n p a
  -> Matrix n p b
  -> Matrix n p c
zip za zb f m1 m2 = izip za zb (\_ _ a b -> f a b) m1 m2

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

{-# INLINE densify #-}
densify
  :: forall (n :: Nat) (p :: Nat).
     (KnownNat n, KnownNat p)
  => Matrix n p Double -> Dense.L n p
densify (Matrix m) = Dense.L . Dense.Dim . Dense.Dim $ M.pack m
