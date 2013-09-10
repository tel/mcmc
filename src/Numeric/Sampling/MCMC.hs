
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

{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE ConstraintKinds #-}

module Numeric.Sampling.MCMC where

import Control.Applicative
import Control.Lens
import Control.Monad
import Control.Monad.Primitive
import Control.Monad.Primitive.Class
import Control.Monad.ST
import Data.Foldable (foldMap)
import Data.PrimRef
import Data.Traversable (Traversable)
import Linear.Affine
import Linear.Vector
import Numeric.AD
import Numeric.AD.Types (AD, Mode)
import Numeric.Log
import Pipes
import System.Random.MWC (Variate, Gen)
import System.Random.MWC.Monad
import System.Random.MWC.Distributions.Monad
import qualified Data.Vector as V

import Linear.V1

-- | Loglik. I'd love to use "Numeric.Log" to tag the output value as
-- log domain, but that won't play nicely with "Numeric.AD".
type LL   f a   = Point f a -> a

-- | A vector field synonym
type VF   f a   = Point f a -> f a

-- | Stepper of a location (a proposal distribution or a final
-- stepping distribution).
type Jump m f a = Point f a -> m (Point f a)

-- | Takes a single, symmetric Metropolis-Hastings step according to
-- the given proposal distribution and goal log-likelihood.

metropolis
  :: (RealFloat a, Precise a, Variate a, Ord a, MonadPrim m) 
     => Jump (Rand m) f a -> LL f a -> Jump (Rand m) f a
metropolis q l x0 =
  do prop <- q x0
     test <- uniformR (0, 1)
     return $ if (test < l prop / l x0) then prop else x0
{-# INLINE metropolis #-}


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

-- | Computes a Hamiltonian step

hamiltonian
  :: (Fractional a, Ord a, Variate a, Traversable f, Additive f, MonadPrim m)
     => Rand m (f a) -- ^ the momentum \"flick\"er. If we had to
                     -- generate this ourselves we'd be far less
                     -- general, both in the constraints to this
                     -- function and to our ability to introduce
                     -- momentum exclusively in certain dimensions,
                     -- for instance.
     -> LL f a -> VF f a -> Int -> a
     -> Jump (Rand m) f a
hamiltonian flick l vf steps eps x0 = do
  momentum <- flick
  let h0 = Hamiltonian x0 momentum eps
      prop = h0 ^?! dropping steps (iterated $ leapfrog vf) . position
  test <- uniformR (0, 1)
  return $ if (test < l prop / l x0) then prop else x0
{-# INLINE hamiltonian #-}

-- | Contract to the compiler stating that 'm' is a bottom layer of a
-- 'MonadPrim' stack.
type IsPrim m = (m ~ BasePrimMonad m, PrimMonad m, MonadPrim m)

-- | Warning, not atomic.
     
modifyPrimRefM
  ::  IsPrim m => PrimRef m a -> (a -> Rand m a) -> Rand m a
modifyPrimRefM var f = modifyPrimRefM' var (liftM (\a -> (a,a)) . f)
{-# INLINE modifyPrimRefM #-}

modifyPrimRefM'
  :: IsPrim m => PrimRef m a -> (a -> Rand m (a, b)) -> Rand m b
modifyPrimRefM' var f = toRand $ \gen -> do
  x <- readPrimRef var
  (x', o) <- runRand (f x) gen
  writePrimRef var x'
  return o
{-# INLINE modifyPrimRefM' #-}

-- | Burn-in then sample `n` points from any jumping algorithm

sample :: IsPrim m => Int -> Int -> Point f a
          -> Jump (Rand m) f a -> Rand m (V.Vector (Point f a))
sample burn n x0 jump = do
  ref <- liftPrim (newPrimRef x0)
  replicateM_ burn $ modifyPrimRefM ref jump
  V.replicateM n $ modifyPrimRefM ref jump
{-# INLINE sample #-}

-- | Reified sampling stream

stream
  :: IsPrim m => Int -> Point f a
     -> Jump (Rand m) f a -> Producer (Point f a) (Rand m) r
stream burn x0 jump = do
  ref <- lift $ liftPrim (newPrimRef x0)
  lift $ replicateM_ burn $ modifyPrimRefM ref jump
  forever $ lift (modifyPrimRefM ref jump) >>= yield


uniQ :: (Num a, Variate a, MonadPrim m) => Jump (Rand m) V1 a
uniQ _ = liftM return $ uniformR (0, 1)

lik :: Floating a => LL V1 a
lik (P (V1 x)) = x
