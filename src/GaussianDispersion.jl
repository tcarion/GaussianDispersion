module GaussianDispersion

using DocStringExtensions
using Parameters: @unpack

export ReleaseParams, DispersionParams, GaussianPlume
export concentration, StabilityClass, Terrain, PasquillGiffordCriteria, pasquill_gifford
export Rural, Urban
# export A, B, C, D, E, F
export Stabilities
export Strong, Moderate, Slight, Cloudy, Clear

const BRIGGS_COEFS = [
    0.22     0.16     0.11     0.08     0.06     0.04    0.2  0.12   0.08     0.06     0.03     0.016    0.32     0.22     0.16     0.11    0.24   0.2   0.14     0.08
    0.0001   0.0001   0.0001   0.0001   0.0001   0.0001  0.0  0.0    0.0002   0.0015   0.0003   0.0003   0.0004   0.0004   0.0004   0.0004  0.001  0.0   0.0003   0.0003
   -0.5     -0.5     -0.5     -0.5     -0.5     -0.5     1.0  1.0   -0.5     -0.5     -1.0     -1.0     -0.5     -0.5     -0.5     -0.5     0.5    1.0  -0.5     -0.5
]

const GRAVITY = 9.80665 

include("meteorology.jl")
include("parametrization.jl")
include("parameters.jl")
include("concentration.jl")

end
