{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
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

dim :: forall (n :: Nat) a. KnownNat n => V n a -> Int
dim _ = fromIntegral n
  where n = natVal (Proxy :: Proxy n)

values :: V n a -> UV.Vector a
values (V v) = V.values v

nnz :: UV.Unbox a => V n a -> Int
nnz = UV.length . values

indices :: V n a -> UV.Vector Int
indices (V v) = V.indices v

singleton :: forall (n :: Nat) a. (KnownNat n, UV.Unbox a) => Int -> a -> V n a
singleton i a = V $ V.Vector n (UV.singleton i) (UV.singleton a)
  where n = fromIntegral $ natVal (Proxy :: Proxy n)

{-# INLINE add #-}
add :: (Num a, UV.Unbox a) => V n a -> V n a -> V n a
add (V a) (V b) = V (a+b)

fromList
  :: forall (n :: Nat) a. (KnownNat n, UV.Unbox a, Num a) => [(Int, a)] -> V n a
fromList xs = V $ n V.|> xs
  where n = fromIntegral $ natVal (Proxy :: Proxy n)

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

{-# INLINE dot #-}
dot :: (Num a, UV.Unbox a) => V n a -> V n a -> a
dot (V u) (V v) = V.dot u v

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
