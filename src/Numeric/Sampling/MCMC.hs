
-- |
-- Module      : Numeric.Sampling.MCMC
-- Copyright   : (c) Joseph Abrahamson 2013
-- License     : MIT
-- 
-- Maintainer  : me@jspha.com
-- Stability   : experimental
-- Portability : non-portable
-- 
-- Symmetric Metropolis sampling.

module Numeric.Sampling.MCMC where

import Linear.Affine
import Control.Monad.Primitive
import Data.Primitive.MutVar
import Control.Monad.ST
import Control.Monad
import System.Random.MWC
import Numeric.Log
import qualified Data.Vector as V

import Linear.V1

-- | Loglik
type LL   f a = Point f a -> Log a

-- | Proposed sample
type Prop m f a = Point f a -> m (Point f a)


-- | Takes a single, symmetric Metropolis-Hastings step according to
-- the given proposal distribution and goal log-likelihood.

metropolis
  :: (RealFloat a, Precise a, Variate a, Ord a, PrimMonad m)
     => Prop m f a -> LL f a -> Gen (PrimState m)
     -> Point f a -> m (Point f a)
metropolis q l gen x0 =
  do prop <- q x0
     test <- uniformR (0, 1) gen
     return $ if (Exp (log test) < l prop / l x0) then prop else x0
{-# INLINE metropolis #-}


-- | Warning, not atomic.
     
modifyMutM :: PrimMonad m => MutVar (PrimState m) a -> (a -> m a) -> m a
modifyMutM var f = do { x <- readMutVar var; x' <- f x; writeMutVar var x'; return x' }
{-# INLINE modifyMutM #-}


-- | Burn-in then sample `n` points using symmetric
-- Metropolis-Hastings steps.

sample
  :: (RealFloat a, Precise a, Variate a, Ord a, PrimMonad m)
     => Prop m f a -> LL f a -> Gen (PrimState m)
     -> Int -> Int -> Point f a -> m (V.Vector (Point f a))
sample q l gen burn n x0 = do
  ref <- newMutVar x0
  replicateM_ burn $ modifyMutM ref (metropolis q l gen)
  V.replicateM n $ modifyMutM ref (metropolis q l gen)
{-# INLINE sample #-}

uniQ :: (Num a, Variate a, PrimMonad m) => Gen (PrimState m) -> Prop m V1 a
uniQ gen _ = liftM return $ uniformR (0, 1) gen

lik :: Floating a => LL V1 a
lik (P (V1 x)) = Exp (log (1-x))
