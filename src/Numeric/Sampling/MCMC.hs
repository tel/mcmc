{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE KindSignatures #-}

module Numeric.Sampling.MCMC where

import Control.Applicative
import Control.Monad
import Control.Monad.Random
import Control.Comonad
import Control.Comonad.Cofree
import Control.Lens
import Linear.Affine
import Linear.Vector
import Numeric.AD
import Numeric.AD.Types (AD, Mode)
import Data.Foldable
import Data.Traversable
import Linear.V2
import System.Random.MWC

import Numeric.Sampling.Util
import Numeric.Sampling.Types

rosen :: RealFloat a => C1Obj V2 a
rosen = naturalC1 $ \(P (V2 x y)) -> (1 - x)^2 + 100 * (y - x ^ 2)^2

naturalC1 :: (Traversable f, RealFloat a) => (forall a . RealFloat a => Obj f a) -> C1Obj f a
naturalC1 p = C1Obj p (unP . grad p)

-- | Compute the initial state space based on a differentiable function.
naturalS
  :: (Num a, Traversable f)
     => (forall s. f (AD s a) -> AD s a)
     -> Point f a -> StateS f a
naturalS f px@(P x) = StateS px (grad f x)

tryReject
  :: (Fractional a, Ord a, Random a, MonadRandom m) =>
     (b -> a) -> b -> b -> m b
tryReject f x0 x1 = do
  coin <- getRandomR (0, 1)
  return (if (coin < f x1 / f x0) then x1 else x0)
{-# INLINE tryReject #-}

-- | Refines an approximate sampler using a metropolis style rejection
-- step. This is only unbiased if the original sampler is symmetric,
-- i.e. @x' == jump x@ is equally likely as @x = jump x'@. This uses
-- direct division in computing the next step, so there may be
-- stability issues unless the base numeric type is compensated.
metropolis
  :: (Fractional a, Ord a, Random a, MonadRandom m)
     => Obj f a -> Sampler m f a -> Sampler m f a
metropolis p jump x0 = jump x0 >>= tryReject p x0
{-# INLINE metropolis #-}

hamiltonian
  :: (MonadRandom m, Fractional a, Ord a, Random a, Traversable f, Additive f)
     => Scalar f a -> Int
     -> Flick m a -> C1Obj f a
     -> Sampler m f a
hamiltonian eps n flick c1 x0 = do
  st <- flickStateS flick x0
  tryReject (c1 ^. fn) x0
    $ st ^?! dropping n (iterated $ leapfrog eps (c1 ^. grd)) . pos
{-# INLINE hamiltonian #-}
