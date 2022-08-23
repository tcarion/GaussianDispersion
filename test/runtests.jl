using GaussianDispersion
using Test

@testset "parameters" begin include("parameters.jl") end
@testset "meteorology" begin include("meteorology.jl") end
@testset "concentration" begin include("concentration.jl") end