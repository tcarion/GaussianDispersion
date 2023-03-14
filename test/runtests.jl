using GaussianDispersion
using Test

const GD = GaussianDispersion

@testset "meteorology" begin include("meteorology.jl") end
@testset "parameters" begin include("parameters.jl") end
@testset "concentration" begin include("concentration.jl") end