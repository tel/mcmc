{-# LANGUAGE RankNTypes #-}

module Numeric.Sampling.DualAveraging where

import Control.Monad
import Control.Monad.Random
import Data.Traversable
import Linear.Vector
import Linear.Affine

import Numeric.Sampling.Util
import Numeric.Sampling.Types

guessEpsilon
  :: (Additive f, Fractional a, MonadRandom m, Traversable f, Random a)
     => Flick m a -> Point f a -> Grad f a -> m (StateS f a)
guessEpsilon flick x0 grad = do
  st <- liftM (leapfrog 1 grad) (flickStateS flick x0)
  let a = 2 * (if )
  return st
