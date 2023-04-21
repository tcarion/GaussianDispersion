using GaussianDispersion
using GaussianDispersion:
    saturation_pressure,
    virtual_temp,
    friction_velocity,
    potential_temp,
    friction_temp,
    obukhov,
    solar_energy_flux,
    net_radiation,
    sensible_heat_flux,
    log_wind_profile,
    log_wind_profile_neutral
using GaussianDispersion: PGStabilityClass, PGStability
using GaussianDispersion: PGSlightlyUnstable, PGNeutral
using GaussianDispersion: AbstractStability
using Test
using Roots

@testset "Pasquill Gifford stability classes" begin
    pgstab = PGStability(:A)
    @test pgstab isa AbstractStability

    @test pasquill_gifford(Moderate(), 5.5) == Set([PGSlightlyUnstable(), PGNeutral()])
end



@testset "scalar atmosphere properties" begin
    cp = 1004.5
    R_gas = 287.05
    ew = saturation_pressure(273.15 + 50.) * 1e-3
    @test ew ≈ 12.3440 atol = 1e-1

    ps = 99922.5703
    t2 = 290.316681
    td2 = 285.705688
    stress = 0.0230113287

    e = saturation_pressure(t2)
    tv = virtual_temp(ps, t2, e)
    u_star = friction_velocity(ps, t2, td2, stress, 287.05)
    @test u_star ≈ 0.138913378 atol = 1e-5
    
    p1 = 99804.1719
    t1 = 291.607758
    θ = potential_temp(291.607758, p1, R_gas, cp)
    @test θ ≈ 291.77113 atol = 1e-4

    q = 6.3768158
    rho = ps / (R_gas * tv)
    θ_star = friction_temp(q, rho, u_star, cp)
    @test θ_star ≈ -0.0383189023 atol = 1e-4

    L = obukhov(ps, t2, td2, t1, p1, u_star, q)
    @test L ≈ -37 atol = 1
end

@testset "heat balance" begin
    
    # Ex 5.4 p.79
    n = 0.
    solar_elev = 60.
    terrain = Rural
    albedo = 0.2
    C_g = 0.1
    B = 0.5
    T_surf = 298.15
    R = solar_energy_flux(solar_elev, n)
    @test R ≈ 827.4 atol = 1e-1
    R_N = net_radiation(R, albedo, n, T_surf)
    @test R_N ≈ 524 atol = 1
    q = sensible_heat_flux(C_g, R_N, B)
    @test q ≈ 157.2 atol = 1
    
end

@testset "solar elevation" begin
    δ = GD.solar_declination(1)
    lat = 52.
    h = 15
    ha = GD.hour_angle(15.84)
    elev = GD.solar_elevation(lat, δ, ha)
    @test elev ≈ 37.73 atol = 1e-2
    @test GD.albedo_elevation(0.1, 30) ≈ 0.129 atol = 1e-2
end

@testset "wind speed profile" begin
    z₀ = 0.1

    # Neutral condition - ex 5.5
    z_m = 10.
    u_m = 7.
    u_star = friction_velocity(u_m, z_m, z₀)
    @test u_star ≈ 0.608 atol=1e-3
    u_90 = log_wind_profile_neutral(u_star, 90., z₀) 
    @test u_90 ≈ 10.3 atol=1e-1

    # Non-Neutral conditions - ex 5.6
    p = 101e3
    T = 300.
    z₀ = 0.5
    q = 200
    u_10 = 5.
    f(u_star) = u_10 - log_wind_profile(u_star, 10., z₀, obukhov(p, T, u_star, q))
    flog(u_star) = log_wind_profile(u_star, 10., z₀, obukhov(p, T, u_star, q))
    u_star = find_zero(f, 0.5)
    @test u_star ≈ 0.711 atol=1e-3
end