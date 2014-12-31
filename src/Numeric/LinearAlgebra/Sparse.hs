{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ViewPatterns #-}

module Numeric.LinearAlgebra.Sparse
    ( CxSparse()
    , Matrix(), nRows, nColumns
    , mul
    , compress, fromTriples, (><)
    , transpose, ctrans, hermitian
    , lin
    , add
    , gaxpy, gaxpy_, mulV
    , hcat, vcat, fromBlocks, fromBlocksDiag
    , kronecker
    , takeDiag, diag, ident
    , zeros
    , module Data.Complex.Enhanced
    ) where

import Control.Applicative
import Control.Monad (void)
import Data.Foldable
import qualified Data.List as List
import Data.MonoTraversable
import Data.Vector.Unboxed (Unbox)
import qualified Data.Vector.Unboxed as U
import Data.Vector.Storable (Vector)
import qualified Data.Vector.Storable as V
import Data.Vector.Storable.Mutable (IOVector)
import qualified Data.Vector.Storable.Mutable as MV
import Foreign.ForeignPtr.Safe (newForeignPtr)
import Foreign.Marshal.Alloc (finalizerFree)
import Foreign.Marshal.Utils (with)
import Foreign.Storable
import GHC.Stack
import Prelude hiding (any, foldl1)
import System.IO.Unsafe (unsafePerformIO)

import Data.Complex.Enhanced
import Data.Matrix.Sparse

mul :: CxSparse a => Matrix a -> Matrix a -> Matrix a
mul = mul_go where
  {-# NOINLINE mul_go #-}
  mul_go _a _b =
    unsafePerformIO $
    withConstCs _a $ \_a ->
    withConstCs _b $ \_b ->
      cs_multiply _a _b >>= fromCs
{-# INLINE mul #-}

compress
  :: (CxSparse a, Unbox a)
  => Int -> Int -> U.Vector (Int, Int, a) -> Matrix a
compress = compress_go where
  {-# NOINLINE compress_go #-}
  compress_go nr nc (U.unzip3 -> (rs, cs, xs)) =
    unsafePerformIO $
    withConstTriples nr nc (V.convert rs) (V.convert cs) (V.convert xs) $ \pcs ->
      cs_compress pcs >>= fromCs
{-# INLINE compress #-}

fromTriples :: CxSparse a => Int -> Int -> [(Int, Int, a)] -> Matrix a
{-# INLINE fromTriples #-}
fromTriples = fromTriples_go where
  {-# NOINLINE fromTriples_go #-}
  fromTriples_go nr nc (unzip3 -> (rs, cs, xs)) =
    unsafePerformIO $
    withConstTriples nr nc (V.fromList rs) (V.fromList cs) (V.fromList xs)
      $ \ptr -> cs_compress ptr >>= fromCs

(><) :: CxSparse a => Int -> Int -> [(Int, Int, a)] -> Matrix a
{-# INLINE (><) #-}
(><) = fromTriples

transpose :: CxSparse a => Matrix a -> Matrix a
transpose = transpose_go where
  {-# NOINLINE transpose_go #-}
  transpose_go mat =
    unsafePerformIO $
    withConstCs mat $ \cs ->
      cs_transpose cs (V.length $ values mat) >>= fromCs
{-# INLINE transpose #-}

ctrans :: (CxSparse a, IsReal a) => Matrix a -> Matrix a
ctrans = omap conj . transpose
{-# INLINE ctrans #-}

hermitian :: (Eq a, IsReal a, CxSparse a) => Matrix a -> Bool
hermitian m = ctrans m == m
{-# INLINE hermitian #-}

lin :: CxSparse a => a -> Matrix a -> a -> Matrix a -> Matrix a
lin = lin_go where
  {-# NOINLINE lin_go #-}
  lin_go _alpha _a _beta _b =
    unsafePerformIO $
    withConstCs _a $ \_a ->
    withConstCs _b $ \_b ->
    with _alpha $ \_alpha ->
    with _beta $ \_beta ->
      cs_add _a _b _alpha _beta >>= fromCs
{-# INLINE lin #-}

add :: (CxSparse a, Num a) => Matrix a -> Matrix a -> Matrix a
add a b = lin 1 a 1 b
{-# INLINE add #-}

gaxpy_ :: (CxSparse a) => Matrix a -> IOVector a -> IOVector a -> IO ()
gaxpy_ _a _x _y =
  withConstCs _a $ \_a ->
  MV.unsafeWith _x $ \_x ->
  MV.unsafeWith _y $ \_y ->
    void $ cs_gaxpy _a _x _y
{-# INLINE gaxpy_ #-}

gaxpy :: CxSparse a => Matrix a -> Vector a -> Vector a -> Vector a
gaxpy = gaxpy_go where
  {-# NOINLINE gaxpy_go #-}
  gaxpy_go a _x _y =
    unsafePerformIO $ do
      _y <- V.thaw _y
      _x <- V.unsafeThaw _x
      gaxpy_ a _x _y
      V.unsafeFreeze _y
{-# INLINE gaxpy #-}

mulV :: (CxSparse a, Num a) => Matrix a -> Vector a -> Vector a
mulV = mulV_go where
  {-# NOINLINE mulV_go #-}
  mulV_go a _x =
    unsafePerformIO $ do
      _x <- V.unsafeThaw _x
      y <- MV.replicate (MV.length _x) 0
      gaxpy_ a _x y
      V.unsafeFreeze y
{-# INLINE mulV #-}

hcat :: Storable a => Matrix a -> Matrix a -> Matrix a
hcat a b
  | nRows a /= nRows b = errorWithStackTrace "row dimension mismatch"
  | otherwise = Matrix
      { nRows = nRows a
      , nColumns = nc
      , columnPointers =
        V.init (columnPointers a) V.++ (V.map (+ nza) $ columnPointers b)
      , rowIndices = rowIndices a V.++ rowIndices b
      , values = values a V.++ values b
      }
  where
    nc = nColumns a + nColumns b
    nza = V.last $ columnPointers a
{-# INLINE hcat #-}

vcat :: (CxSparse a, Storable a) => Matrix a -> Matrix a -> Matrix a
vcat a b
  | nColumns a /= nColumns b = errorWithStackTrace "column dimension mismatch"
  | otherwise = transpose $ hcat (transpose a) (transpose b)
{-# INLINE vcat #-}

fromBlocks :: (CxSparse a, Storable a) => [[Matrix a]] -> Matrix a
fromBlocks = foldl1 vcat . map (foldl1 hcat)
{-# INLINE fromBlocks #-}

fromBlocksDiag :: (CxSparse a, Storable a) => [[Matrix a]] -> Matrix a
fromBlocksDiag = foldl1 vcat . zipWith joinAt [0..] . List.transpose where
  {-# INLINE joinAt #-}
  joinAt = \n mats ->
    case splitAt n mats of
     ([], ls) -> foldl1 hcat ls
     (rs, []) -> foldl1 hcat rs
     (rs, ls) -> foldl1 hcat ls `hcat` foldl1 hcat rs
{-# INLINE fromBlocksDiag #-}

kronecker :: CxSparse a => Matrix a -> Matrix a -> Matrix a
kronecker = kronecker_go where
  {-# NOINLINE kronecker_go #-}
  kronecker_go _a _b =
    unsafePerformIO $
    withConstCs _a $ \_a ->
    withConstCs _b $ \_b ->
      cs_kron _a _b >>= fromCs
{-# INLINE kronecker #-}

takeDiag :: CxSparse a => Matrix a -> Vector a
takeDiag = takeDiag_go where
  {-# NOINLINE takeDiag_go #-}
  takeDiag_go _a@Matrix{..} =
    unsafePerformIO $
    withConstCs _a $ \_a ->
      V.unsafeFromForeignPtr0
      <$> (cs_diag _a >>= newForeignPtr finalizerFree)
      <*> pure (min nRows nColumns)
{-# INLINE takeDiag #-}

diag :: Storable a => Vector a -> Matrix a
diag values = Matrix{..}
  where
    nColumns = V.length values
    nRows = nColumns
    columnPointers = V.iterateN (nRows + 1) (+1) 0
    rowIndices = V.iterateN nColumns (+1) 0
{-# INLINE diag #-}

ident :: (Num a, Storable a) => Int -> Matrix a
ident n = diag $ V.replicate n 1
{-# INLINE ident #-}

zeros :: (Num a, Storable a) => Int -> Int -> Matrix a
zeros nRows nColumns = Matrix{..}
  where
    columnPointers = V.replicate (nColumns + 1) 0
    rowIndices = V.empty
    values = V.empty
{-# INLINE zeros #-}
