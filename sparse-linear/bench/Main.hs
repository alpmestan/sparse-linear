{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Main where

import qualified Data.Matrix.Sparse.Static    as SM
import qualified Data.Vector.Sparse.Static    as SV
import qualified Numeric.LinearAlgebra.Static as D

import Criterion.Main
import GHC.TypeLits
import Data.Proxy
import Data.Reflection
import System.Random
import Control.Monad (replicateM)
import Data.Traversable (forM)
import Debug.Trace

main :: IO ()
main = defaultMain
  [ bgroup "mXv"
    [ randMat n $ \mat -> bgroup ("n = " ++ show n)
      [ bench "spmv" $ whnf f_spmv mat
      , bench "spmspv1" $ whnf f_spmspv1 mat
      ]
    | n <- [ 10, 100, 1000, 2000, 4000, 10000 ]
    ]
  ]

randMat :: Int -> (forall n. KnownNat n => SM.Matrix n n Double -> Benchmark) -> Benchmark
randMat n f = env randMat' $ \xs ->
  reifyNat (fromIntegral n) $
    \(_pn :: Proxy n) -> f (SM.fromList xs :: SM.Matrix n n Double)

  where randMat' = fmap concat $
          forM [0..(n-1)] $ \j ->
            map (\(i, a) -> (i, j, a)) <$> randCol
        randCol  = zip [0..(k-1)] <$> replicateM k randomIO
        k = max 5 $ (n-1) `div` 100

v_d :: forall (n :: Nat). KnownNat n => D.R n
v_d = D.vector $ 1 : replicate (n-1) 0

  where n = fromIntegral $ natVal (Proxy :: Proxy n)

v_s :: forall (n :: Nat). KnownNat n => SV.V n Double
v_s = SV.singleton 0 1.0

f_spmv :: KnownNat n => SM.Matrix n n Double -> D.R n
f_spmv m = iterate (SM.mulV m) v_d !! 3

f_spmspv1 :: forall n. KnownNat n => SM.Matrix n n Double -> SM.Matrix n 1 Double
f_spmspv1 m = iterate (SM.mul m) (SM.asColumn v_s) !! 3
