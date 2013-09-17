{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE KindSignatures #-}

module Numeric.Sampling.MCMC2 where

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
import Numeric.AD.Internal.Classes (Lifted)
import Data.Foldable
import Data.Traversable
import Linear.V2
import System.Random.MWC

type Scalar    (f :: * -> *) a =                         a
type Obj       (f :: * -> *) a = Point f a ->            a
type Grad      (f :: * -> *) a = Point f a ->          f a
type Sampler m (f :: * -> *) a = Point f a -> m (Point f a)

data StateS f a = StateS { _pos :: Point f a, _mom :: f a }
makeLenses ''StateS

data C1Obj f a = C1Obj { _fn :: Obj f a, _grd :: Grad f a }
makeLenses ''C1Obj

class (RealFrac a, Floating a) => Smooth a
instance (Lifted f, Smooth a) => Smooth (AD f a)

unP :: Point f x -> f x
unP (P x) = x
{-# INLINE unP #-}

rosen :: Smooth a => C1Obj V2 a
rosen = naturalC1 $ \(P (V2 x y)) -> (1 - x)^2 + 100 * (y - x ^ 2)^2

naturalC1 :: (Traversable f, Smooth a) => (forall a . Smooth a => Obj f a) -> C1Obj f a
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

hamiltonian
  :: (MonadRandom m, Smooth a, Ord a, Random a, Traversable f, Additive f)
     => Scalar f a -> Int
     -> (Point f a -> m (StateS f a)) -> C1Obj f a
     -> Sampler m f a
hamiltonian eps n mkStM c1 x0 = do
  st <- mkStM x0
  tryReject (c1 ^. fn) x0
    $ st ^?! dropping n (iterated $ leapfrog eps (c1 ^. grd)) . pos
{-# INLINE hamiltonian #-}

-- | Not particularly great random momentum "flicking"
randomMomentum
  :: (Traversable f, MonadRandom m, Random a)
     => Point f a -> m (StateS f a)
randomMomentum x = liftM (StateS x . unP) (Data.Traversable.mapM (const getRandom) x)
