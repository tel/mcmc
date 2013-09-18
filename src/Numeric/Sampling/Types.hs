{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE KindSignatures #-}

module Numeric.Sampling.Types where

import Control.Lens
import Linear.Affine
import Numeric.AD.Internal.Classes (Lifted)
import Numeric.AD.Types (AD, Mode)

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