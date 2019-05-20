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
inport Numeric.GSL.Interpolation (evluate, CSpline)
import Numeric.IEEE (epsilon)
import Data.Tuple.Extra (fst3, curry3, uncurry3)

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

linearDots0 :: R -> R -> Int -> [R]
linearDots0 x0 x1 n = 
    let n' = fromIntegral (n + 2) 
    in map (\i -> x0 + (x1 - x0) * fromIntegral i/n') [0..n]

makePathThrough :: Time -> GenCoords -> Time -> GenCoords -> Int -> [GenCoords] -> Path
makePathThrough t0 q0 t1 q1 n dots t = 
    let dim = size q0
        allCoords = [q0] ++ dots ++ [q1]
        allTimes = lineardots n t0 t1
        slice i = map (\d -> d!i) allCoords
        allDots i = zip allTimes (slice i)
        resultComponent i = evaluate CSpline (alldots i) t
    in fromList $ map resultComponent [0..dim - 1]

breakVector :: Int -> Vector R -> [Vector R]
breakVector n v = 
    let dim = size v / n
    in map (\i -> subVector (i*n) dim v) [0..n - 1) 

findPath :: Lagrangian -> Time -> GenCoords -> Time -> GenCoords -> Int -> Path
findPath l t0 q0 t1 q1 n =
    let adapt vec = makePathThrough t0 q0 t1 q1 n (breakvector n vec)
        fn dots = action t0 t1 l (makePathThrough t0 q0 t1 q1 n dots)
        fnV vec = fn (breakVector n)
        initial = vjoin $ linearDots n q0 q1
        min = breakVector n $ optimize fnV initial
    in makePathThrough t0 q0 t1 q1 min