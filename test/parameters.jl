using GaussianDispersion
using GaussianDispersion: Rural, Urban
using GaussianDispersion: HeatBalanceParams, pasquill_gifford
using GaussianDispersion: Aermod, Calpuff
using Test

@testset "dispersion params" begin
    terrain = Rural()
    stability = GaussianDispersion.A

    @test DispersionParams(0., terrain, stability) == DispersionParams(0., 0.)
end

@test pasquill_gifford(Moderate, 5.5) == Stabilities(:C, :D)
@test GaussianDispersion.briggs_dispersion(Urban, GaussianDispersion.D) == (
    (0.16, 0.0004, -0.5),
    (0.14, 0.0003, -0.5)
    )

@test HeatBalanceParams(Aermod) isa HeatBalanceParams
@test HeatBalanceParams(C_g = 0.3) isa HeatBalanceParams