
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

{-# LANGUAGE RankNTypes, TemplateHaskell #-}

module Numeric.Sampling.MCMC where

import Control.Applicative
import Control.Lens
import Control.Monad
import Control.Monad.Primitive
import Control.Monad.ST
import Data.Primitive.MutVar
import Data.Traversable (Traversable)
import Linear.Affine
import Linear.Vector
import Numeric.AD
import Numeric.AD.Types (AD, Mode)
import Numeric.Log
import System.Random.MWC
import qualified Data.Vector as V

import Linear.V1

-- | Loglik. I'd love to use "Numeric.Log" to tag the output value as
-- log domain, but that won't play nicely with "Numeric.AD".
type LL   f a = Point f a -> a

-- | A vector field synonym
type VF   f a = Point f a -> f a

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
     return $ if (test < l prop / l x0) then prop else x0
{-# INLINE metropolis #-}


hamiltonian
  :: (RealFloat a, Precise a, Variate a, Ord a, PrimMonad m)
     => Prop m f a -> LL f a -> Gen (PrimState m)
     -> Point f a -> m (Point f a)
hamiltonian = undefined

data Hamiltonian f a =
  Hamiltonian { _position :: Point f a
              , _momentum :: f a
              , _epsilon  :: a
              }
makeLenses ''Hamiltonian

-- | Declares that a point should be interpreted as a vector.

isVec :: Point f a -> f a
isVec (P a) = a

(^+^~) :: (Num a, Additive f) => Setting (->) s t (f a) (f a) -> f a -> s -> t
l ^+^~ a = over l (^+^ a)
{-# INLINE (^+^~) #-}

(.+^~) :: (Num a, Affine f) => Setting (->) s t (f a) (f a) -> Diff f a -> s -> t
l .+^~ a = over l (.+^ a)
{-# INLINE (.+^~) #-}

(^-^~) :: (Num a, Additive f) => Setting (->) s t (f a) (f a) -> f a -> s -> t
l ^-^~ a = over l (^-^ a)
{-# INLINE (^-^~) #-}


-- | 'grad' with a slightly more informative type signature.

gradP
  :: (Num a, Traversable f)
     => (forall s. Mode s => Point f (AD s a) -> AD s a)
     -> Point f a -> f a
gradP ll = isVec . grad ll
{-# INLINE gradP #-}

-- | Perform a single leapfrog integration step

leapfrog
  :: (Traversable f, Additive f, Fractional a)
     => VF f a -> Hamiltonian f a -> Hamiltonian f a
leapfrog g = leap g . frog . leap g
  where
    half :: Fractional a => a -> a
    half x = x/2
    {-# INLINE half #-}

    leap :: (Additive f, Fractional a) => VF f a -> Hamiltonian f a -> Hamiltonian f a
    leap gradAt h = h & momentum ^+^~ (h ^. epsilon . to half . to (*^ perturbation))
      where perturbation = h ^. position . to gradAt
    {-# INLINE leap #-}

    frog :: (Additive f, Fractional a) => Hamiltonian f a -> Hamiltonian f a
    frog h = h & position .+^~ (eps *^ mom) where
      mom = h ^. momentum
      eps = h ^. epsilon . to half
    {-# INLINE frog #-}

{-# INLINE leapfrog #-}



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
lik (P (V1 x)) = (1-x)
