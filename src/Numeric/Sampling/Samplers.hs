
module Numeric.Sampling.Samplers where

import Pipes
import qualified Pipes.Prelude as P
import Linear.Affine

import Numeric.Sampling.Types

-- | This is nothing more than a type specialized "iterateM"
-- 
--     iterateM :: (a -> m a) -> a -> Producer a m r
-- 
streamS :: Monad m => Sampler m f a -> Point f a -> Producer (Point f a) m r
streamS s x0 = go x0 where go x = yield x >> lift (s x) >>= go
