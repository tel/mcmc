{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RankNTypes #-}

module Numeric.Sampling.Util where

import Control.Lens
import Control.Monad
import Control.Monad.Random
import Data.Traversable
import Linear.Affine
import Linear.Vector

import Numeric.Sampling.Types

unP :: Point f x -> f x
unP (P x) = x
{-# INLINE unP #-}

(^+^~) :: (Num a, Additive f) => Setting (->) s t (f a) (f a) -> f a -> s -> t
l ^+^~ a = over l (^+^ a)
{-# INLINE (^+^~) #-}

(.+^~) :: (Num a, Affine f) => Setting (->) s t (f a) (f a) -> Diff f a -> s -> t
l .+^~ a = over l (.+^ a)
{-# INLINE (.+^~) #-}


-- | A Hamiltonian leapfrog integration step.
leapfrog :: (Additive f, Fractional a) => Scalar f a -> Grad f a -> StateS f a -> StateS f a
leapfrog eps field = leap . frog . leap where
  leap h = h & mom ^+^~ (eps/2 *^ field (h ^. pos))
  {-# INLINE leap #-}
  frog h = h & pos .+^~ (eps/2 *^       (h ^. mom))
  {-# INLINE frog #-}
{-# INLINE leapfrog #-}

flickStateS
  :: (Traversable f, MonadRandom m, Random a)
     => Flick m a -> Point f a -> m (StateS f a)
flickStateS flick x = liftM (StateS x . unP) (flick x)

-- | Not particularly great random momentum "flicking"
uniformMomentumStateS
  :: (Traversable f, MonadRandom m, Random a)
     => Point f a -> m (StateS f a)
uniformMomentumStateS x = liftM (StateS x . unP) (uniformFlick x)

uniformFlick :: (MonadRandom m, Random a) => Flick m a
uniformFlick = Data.Traversable.mapM (const getRandom)
