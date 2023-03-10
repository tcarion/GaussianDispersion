using GaussianDispersion
using GaussianDispersion: Rural, Urban
using GaussianDispersion: HeatBalanceParams, pasquill_gifford
using GaussianDispersion: Aermod, Calpuff
using GaussianDispersion: PGStability
using GaussianDispersion: AbstractDispersionFunctions
using GaussianDispersion: BriggsFunctions, sigma_y, sigma_z
using GaussianDispersion: DispersionCoefficients
using Test

@testset "Disperions parametrization" begin
    struct Foo <: AbstractDispersionFunctions end

    f = Foo()
    @test f(3.) == (3., 3.)
end

@testset "Briggs parametrization" begin
    briggsA = BriggsFunctions(Rural(), PGStability(:A))

    briggsB = BriggsFunctions(Urban(), PGStability(:B))
    # julia> @btime BriggsFunctions(Rural(), PGStability(:A))(5.)
    # 544.300 ns (16 allocations: 496 bytes)
    # (1.09972510308205, 1.0)

    # julia> @btime (0.22 * 5. * (1 + -0.0001*5.)^-0.5, 0.2 * 5. * (1 + 0.0*5.)^1.)
    # 1.495 ns (0 allocations: 0 bytes)
    # (1.1002751031679876, 1.0)

    # julia> @btime briggsA(5.)
    # 204.316 ns (5 allocations: 96 bytes)
    # (1.09972510308205, 1.0)

    syA = sigma_y(briggsA)
    syB = sigma_y(briggsB)

    @test syA(4.) < syB(4.) 

    coefs = DispersionCoefficients(briggsA, 4.)

    @test coefs.y == syA(4.) && coefs.z == sigma_z(briggsA)(4.)
end

@testset "heat balance params" begin
    @test HeatBalanceParams(Aermod) isa HeatBalanceParams
    @test HeatBalanceParams(C_g = 0.3) isa HeatBalanceParams

    @test HeatBalanceParams(Calpuff, Urban).albedo == 0.18
end