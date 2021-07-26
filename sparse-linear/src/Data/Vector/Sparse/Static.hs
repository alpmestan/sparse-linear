{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Data.Vector.Sparse.Static where

import qualified Data.List           as L
import qualified Data.Vector.Sparse  as V
import qualified Data.Vector.Unboxed as UV
import GHC.TypeLits
import qualified Data.Vector.Storable as VS
import Numeric.LinearAlgebra.Static (R, withVector)
import Unsafe.Coerce (unsafeCoerce)
import Data.Proxy

newtype V (n :: Nat) a = V (V.Vector UV.Vector a)
  deriving (Eq, Show)

add :: (Num a, UV.Unbox a) => V n a -> V n a -> V n a
add (V a) (V b) = V (V.glin 0 (+) a (+) b)

sum :: forall (n :: Nat) a. (Num a, UV.Unbox a, KnownNat n) => [V n a] -> V n a
sum xs = L.foldl' add z xs
  where z = V (V.Vector n UV.empty UV.empty) :: V n a
        n = fromIntegral $ natVal (Proxy :: Proxy n)

map :: (UV.Unbox a, UV.Unbox b) => (a -> b) -> V n a -> V n b
map f (V a) = V (V.cmap f a)

(!) :: (UV.Unbox a, Num a) => V n a -> Int -> a
(!) (V v) i = v V.! i

sumV :: (Num a, UV.Unbox a) => V n a -> a
sumV (V v) = UV.sum (V.values v)

toDense :: V n Double -> R n
toDense (V v) = withVector vec unsafeCoerce
  where go (i, i_nnz)
          | i == n = Nothing
          | i_nnz == n_nnz = Just (0, (i+1, i_nnz))
          | i == (ids UV.! i_nnz) = Just (vals UV.! i_nnz, (i+1, i_nnz+1))
          | otherwise = Just (0, (i+1, i_nnz))
        n = V.length v
        n_nnz = UV.length ids
        ids = V.indices v
        vals = V.values v
        vec = VS.unfoldr go (0, 0)
