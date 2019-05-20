module Lagrange
    ( plotPath
    , minmax
    , eps
    ) where

import Control.Exception (assert)
import Data.List
import qualified Graphics.Gnuplot.Simple as GPS
import Numeric.LinearAlgebra.Data
import Numeric.GSL.Differentiation (derivCentral)
import Numeric.GSL.Integration (integrateCQUAD)
import Numeric.GSL.Interpolation (evaluate, InterpolationMethod(..))
import Numeric.GSL.Minimization (minimize, MinimizeMethod(..))
import Numeric.IEEE (epsilon)
import Data.Tuple.Extra (fst3, curry3, uncurry3)
import Data.List.Extra (chunksOf)

type Time = R
type GenCoords = Vector R
type Path = Time -> GenCoords
type LocalTuple = (Time, GenCoords, GenCoords)
type Lagrangian = LocalTuple -> R

minmax :: Ord a => [a] -> (a, a)
minmax (x:xs) = foldl' (\(min', max') x -> (min min' x, max max' x)) (x, x) xs

plotPath :: [(Double, Double)] -> IO ()
plotPath dots = 
            GPS.plotList [] dots

basisVector :: Int -> Int -> Vector R
basisVector n i = 
    let zeros = repeat 0
    in fromList $ concat [take i zeros, [1], take (n - i -1) zeros]

eps :: R
eps = sqrt epsilon

differential :: Path -> R -> Vector R
differential p t =
    let ptP = p (t + eps)
        ptM = p (t - eps)
    in (ptP - ptM) / (scalar $ 2 * eps)

integrate :: (R -> R) -> R -> R -> R
integrate f a b = fst3 $ integrateCQUAD eps 1000 f a b

gamma :: Path -> R -> LocalTuple
gamma p t = 
    let c = p t
        v = differential p t
    in (t, c, v)

action :: R -> R -> Lagrangian -> Path -> R
action t0 t1 l p = integrate (\t -> l (gamma p t)) t0 t1

harmonicOscillator :: R -> R -> Lagrangian
harmonicOscillator m k (t, x, v) = 
    0.5 * m * (v!0)**2 - 0.5 * k * (x!0)**2

linearDots :: Int -> GenCoords -> GenCoords -> [GenCoords]
linearDots n q0 q1 = error "todo"

equalSplit :: Int -> R -> R -> [R]
equalSplit n x0 x1 = 
    let n' = fromIntegral (n + 2) 
    in map (\i -> x0 + (x1 - x0) * fromIntegral i/n') [0..n]

-- TODO: this re-interpolates for each t
makePathThrough :: Time -> GenCoords -> Time -> GenCoords -> Int -> [GenCoords] -> Path
makePathThrough t0 q0 t1 q1 n dots t = 
    let dim = size q0
        allCoords = [q0] ++ dots ++ [q1]
        allTimes = equalSplit n t0 t1
        slice i = map (\d -> d!i) allCoords
        allDots i = zip allTimes (slice i)
        resultComponent i = evaluate CSpline (allDots i) t
    in fromList $ map resultComponent [0..dim - 1]

findPath :: Lagrangian -> Time -> GenCoords -> Time -> GenCoords -> Int -> Path
findPath l t0 q0 t1 q1 n =
    let dim = size q0
        list2Vecs = map fromList . chunksOf dim
        pathThrough = makePathThrough t0 q0 t1 q1 n . list2Vecs
        fn = action t0 t1 l . pathThrough
        initial = concat $ map toList $ linearDots n q0 q1
        steps = take (dim*n) $ repeat 0.1
        minVec = fst $ minimize NMSimplex eps 1000 steps fn initial
    in pathThrough minVec