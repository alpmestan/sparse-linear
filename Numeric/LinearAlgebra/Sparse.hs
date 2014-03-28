{-# LANGUAGE TypeFamilies #-}

module Numeric.LinearAlgebra.Sparse where

import Control.Applicative ((<$>))
import Control.Monad (when)
import Control.Monad.Primitive (PrimMonad(..))
import Control.Monad.ST (runST)
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Unboxed.Mutable as MU
import Data.Vector.Algorithms.Intro (sortBy)
import Data.Vector.Algorithms.Search (Comparison, binarySearchP)

data CS a = CS !Int !(U.Vector Int) !(U.Vector (Int, a))
data Coord a = Coord !Int !Int !(U.Vector (Int, Int, a))

data CSR = CSR
data CSC = CSC
data COO = COO

type family Matrix fmt a
type instance Matrix CSR a = CS a
type instance Matrix CSC a = CS a
type instance Matrix COO a = Coord a

class Format f where
    coo :: Matrix f a -> Matrix COO a
    csr :: Matrix f a -> Matrix CSR a
    csc :: Matrix f a -> Matrix CSC a

rowOrdering, colOrdering :: Comparison (Int, Int, a)
rowOrdering (ra, ca, _) (rb, cb, _) =
    case compare ra rb of
      EQ -> compare ca cb
      x -> x
colOrdering (ra, ca, _) (rb, cb, _) =
    case compare ca cb of
      EQ -> compare ra rb
      x -> x

instance Format COO where
    coo = id

    csr (Coord nRows nCols vals) = runST $ do
      valsM <- U.thaw vals
      sortBy rowOrdering valsM
      rows <- U.generateM nRows
        $ \i -> binarySearchP (\(r, _, _) -> r == i) valsM
      vals' <- U.map (\(_, c, x) -> (c, x)) <$> U.freeze valsM
      return $! CS nCols rows vals'

    csc (Coord nRows nCols vals) = runST $ do
      valsM <- U.thaw vals
      sortBy colOrdering valsM
      cols <- U.generateM nRows
        $ \i -> binarySearchP (\(_, c, _) -> c == i) valsM
      vals' <- U.map (\(r, _, x) -> (r, x)) <$> U.freeze valsM
      return $! CS nRows cols vals'

instance Format CSR where
    coo (CS nCols rows vals) =
        Coord nRows nCols $ U.create $ do
          let spans = U.postscanr' (\start (end, _) -> (start, end)) (len, 0) rows
          vals' <- U.unsafeThaw $ U.map (\(c, x) -> (-1, c, x)) vals
          U.forM_ (U.indexed spans) $ \(r, (start, end)) ->
            forFromUpto start end $ \i -> do
              (_, c, x) <- MU.read vals' i
              MU.write vals' i (r, c, x)
          return vals'
      where
        nRows = U.length rows
        len = U.length vals

    csr = id
    csc = csc . coo

forFromUpto :: Monad m => Int -> Int -> (Int -> m a) -> m a
forFromUpto start end act = go start
  where
    go i | i < end = act i >> go (succ i)
         | otherwise = return ()

instance Format CSC where
    coo (CS nRows cols vals) =
        Coord nRows nCols $ U.create $ do
          let spans = U.postscanr' (\start (end, _) -> (start, end)) (len, 0) cols
          vals' <- U.unsafeThaw $ U.map (\(r, x) -> (r, -1, x)) vals
          U.forM_ (U.indexed spans) $ \(c, (start, end)) ->
            forFromUpto start end $ \i -> do
              (r, _, x) <- MU.read vals' i
              MU.write vals' i (r, c, x)
          return vals'
      where
        nCols = U.length cols
        len = U.length vals

    csr = csr . coo
    csc = id

mm :: (Format f, Format g) => Matrix f a -> Matrix g a -> Matrix COO a
mv :: (Format f, U.Vector v a) => Matrix f a -> v a -> v a

mvM :: (Format f, PrimMonad m, MU.MVector v a)
    => Matrix f a -> v (PrimState m) a -> m (v (PrimState m) a)
