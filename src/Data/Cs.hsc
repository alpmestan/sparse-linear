{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE RecordWildCards #-}

module Data.Cs (Cs(..), CInt(..), Ptr, CxSparse(..)) where

import Control.Applicative
import Data.Vector.Unboxed (Unbox)
import Foreign.C.Types (CInt(..))
import Foreign.Ptr (Ptr)
import Foreign.Storable

import Data.Complex.Enhanced

#include "cs.h"

-- | Translation of the sparse matrix structs used by CxSparse. May
-- contain a sparse matrix in triplet format or compressed sparse column
-- format. In triplet format, the arrays 'p', 'i', and 'x' all have size
-- 'nzmax'; the actual number of entries is given by 'nz'. In compressed
-- format, 'p' is smaller (size 'n + 1') and 'nz == -1'; the actual
-- number of entries is given by the last element of 'p'.
data Cs a = Cs
  { nzmax :: !CInt -- ^ maximum number of entries
  , m :: !CInt -- ^ number of rows
  , n :: !CInt -- ^ number of columns
  , p :: !(Ptr CInt) -- ^ column indices (size nzmax) or column pointers (size n+1)
  , i :: !(Ptr CInt) -- ^ row indices (size nzmax)
  , x :: !(Ptr a) -- ^ numerical values (size nzmax)
  , nz :: !CInt -- ^ number of entries (triplet format) or -1 (compressed format)
  }

instance Storable (Cs (Complex Double)) where
  sizeOf _ = #size cs_ci
  alignment _ = #size cs_ci
  peek ptr = Cs
    <$> (#peek cs_ci, nzmax) ptr
    <*> (#peek cs_ci, m) ptr
    <*> (#peek cs_ci, n) ptr
    <*> (#peek cs_ci, p) ptr
    <*> (#peek cs_ci, i) ptr
    <*> (#peek cs_ci, x) ptr
    <*> (#peek cs_ci, nz) ptr
  poke ptr Cs{..} = do
    (#poke cs_ci, nzmax) ptr nzmax
    (#poke cs_ci, m) ptr m
    (#poke cs_ci, n) ptr n
    (#poke cs_ci, p) ptr p
    (#poke cs_ci, i) ptr i
    (#poke cs_ci, x) ptr x
    (#poke cs_ci, nz) ptr nz
  {-# INLINE sizeOf #-}
  {-# INLINE alignment #-}
  {-# INLINE peek #-}
  {-# INLINE poke #-}

instance Storable (Cs Double) where
  sizeOf _ = #size cs_di
  alignment _ = #size cs_di
  peek ptr = Cs
    <$> (#peek cs_di, nzmax) ptr
    <*> (#peek cs_di, m) ptr
    <*> (#peek cs_di, n) ptr
    <*> (#peek cs_di, p) ptr
    <*> (#peek cs_di, i) ptr
    <*> (#peek cs_di, x) ptr
    <*> (#peek cs_di, nz) ptr
  poke ptr Cs{..} = do
    (#poke cs_di, nzmax) ptr nzmax
    (#poke cs_di, m) ptr m
    (#poke cs_di, n) ptr n
    (#poke cs_di, p) ptr p
    (#poke cs_di, i) ptr i
    (#poke cs_di, x) ptr x
    (#poke cs_di, nz) ptr nz
  {-# INLINE sizeOf #-}
  {-# INLINE alignment #-}
  {-# INLINE peek #-}
  {-# INLINE poke #-}

type CsMultiply a = Ptr (Cs a) -> Ptr (Cs a) -> IO (Ptr (Cs a))

class (Num a, Storable a, Storable (Cs a), Unbox a) => CxSparse a where
  cs_multiply :: CsMultiply a

foreign import ccall "cs.h cs_ci_multiply"
  cs_ci_multiply :: CsMultiply (Complex Double)

instance CxSparse (Complex Double) where
  {-# INLINE cs_multiply #-}
  cs_multiply = cs_ci_multiply

foreign import ccall "cs.h cs_di_multiply"
  cs_di_multiply :: CsMultiply Double

instance CxSparse Double where
  {-# INLINE cs_multiply #-}
  cs_multiply = cs_di_multiply