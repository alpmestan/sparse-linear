{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE RecordWildCards #-}

module Data.Matrix.Sparse.Foreign
       ( fromCs, withConstCs
       ) where

import Control.Applicative
import qualified Data.Vector.Storable as V
import Foreign.ForeignPtr.Safe (newForeignPtr)
import Foreign.Marshal.Alloc (free, finalizerFree)
import Foreign.Marshal.Utils (with)
import Foreign.Ptr
import Foreign.Storable
import GHC.Stack (errorWithStackTrace)

import Data.Cs
import Data.Matrix.Sparse.Compress
import Data.Matrix.Sparse.Type

fromCs :: CxSparse a => Bool -> Ptr (Cs a) -> IO (Matrix Col a)
fromCs doDedup _ptr
  | _ptr == nullPtr = errorWithStackTrace "fromCs: null pointer"
  | otherwise = do
      Cs{..} <- peek _ptr
      let nr = fromIntegral m
          nc = fromIntegral n
          nzmax_ = fromIntegral nzmax
          mkVector ptr len =
            V.unsafeFromForeignPtr0
            <$> newForeignPtr finalizerFree ptr
            <*> pure len
      _rows <- mkVector i nzmax_
      _vals <- mkVector x nzmax_
      _cols <- mkVector p (if nz < 0 then nc + 1 else nzmax_)
      mat <- if nz < 0
             then if doDedup
                  then return $ compress nr nc _rows (decompress _cols) _vals
                  else return Matrix
                    { dimM = nc
                    , dimN = nr
                    , pointers = _cols
                    , indices = _rows
                    , values = _vals
                    }
             else return $ compress nr nc _rows _cols _vals
      free _ptr
      return mat

withConstCs :: CxSparse a => Matrix Col a -> (Ptr (Cs a) -> IO b) -> IO b
withConstCs Matrix{..} act = do
    let nzmax = fromIntegral $ V.length values
        m = fromIntegral dimN
        n = fromIntegral dimM
        nz = -1
    V.unsafeWith pointers $ \p ->
      V.unsafeWith indices $ \i ->
      V.unsafeWith values $ \x ->
      with Cs{..} act