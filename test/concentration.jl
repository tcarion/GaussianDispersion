using GaussianDispersion
using GaussianDispersion.Accessors: @set
# Example from reference book "Air Dispersion Modeling - Foundations and Applications, Alex De Visscher, p.24-25"
@testset "concentration" begin
    # Example 2.1
    u = 7.
    h = 75 + 15
    stabs = pasquill_gifford(Moderate(), u)
    relpar = ReleaseParams(h = h, Q = 100)
    meteopar = MeteoParams(wind = u, stability = collect(stabs)[1])
    plume = GaussianPlume(release = relpar, meteo = meteopar)
    cground1 = plume(1500, 0, 0) * 1e6
    ccross = plume(1500, 100, 0) * 1e6

    @test cground1 ≈ 160.3 atol=1e-1
    @test ccross ≈ 107.5 atol=1e-1

    h = 30
    relpar = ReleaseParams(h = h, Q = 5)

    meteopar = @set meteopar.stability = PGStability(:D)
    meteopar = @set meteopar.wind = 2.
    plume = GaussianPlume(release = relpar, meteo = meteopar)
    
    plume = @set plume.reflection = true

    cground2 = plume(1000, 0, 0)
    @test cground2 ≈ 201.139895390322 * 1e-6

    # julia> @btime result = [plume(x,y,z) for x in range(0, 2000, 200), y in range(-300, 300, 100), z in range(0, 150, 150)]
    # 2.912 s (72000005 allocations: 1.54 GiB)
end

# Example from reference book "Air Dispersion Modeling - Foundations and Applications, Alex De Visscher, p.31"
@testset "Plume rise" begin
    Q = 20 # [m^3/s]
    ρₐ = 1.17
    ρₛ = 0.935
    rₛ = 1.
    wₛ = Q / (π*rₛ^2)
    u = 3

    Fb = GaussianDispersion.buoyancy_flux(ρₛ, ρₐ, rₛ, wₛ)
    @test Fb ≈ 12.539 atol=1e-3
    
    Δh = GaussianDispersion.plume_rise(Fb, 1000., u)
    @test Δh ≈ 47.589 atol=1e-3
end

@testset "concentration - multiple classes" begin
    meteoparams = MeteoParams(
        stability = [PGStability(:A), PGStability(:B)]
    )
    plumeboth = GaussianPlume(meteo = meteoparams)
    cboth= plumeboth(1000, 0, 0)

    plumeA = GaussianPlume(meteo = MeteoParams(stability = PGStability(:A)))
    ca= plumeA(1000, 0, 0)

    plumeB = GaussianPlume(meteo = MeteoParams(stability = PGStability(:B)))
    cb= plumeB(1000, 0, 0)

    @test ca < cboth < cb

    # ! Taking the mean of dispersion functions seems to introduce a huge overhead (x4 computing time !!). I don't know it this is normal...
    # julia> @btime result = [plumeboth(x,y,z) for x in range(0, 2000, 200), y in range(-300, 300, 100), z in range(0, 150, 150)]
    # 9.694 s (153000005 allocations: 4.27 GiB)

    # julia> @btime result = [plumeA(x,y,z) for x in range(0, 2000, 200), y in range(-300, 300, 100), z in range(0, 150, 150)];
    # 2.562 s (75000005 allocations: 1.59 GiB)

    # ! This is weird since:
    # julia> ps = rand(3)
    # julia> @btime (_briggs_function(ps...)(5.), _briggs_function(ps...)(5.))
    # 445.333 ns (11 allocations: 224 bytes)
    # (2.78303613775333, 2.78303613775333)

    # julia> @btime mean([_briggs_function(ps...)(5.), _briggs_function(ps...)(5.)])
    # 549.160 ns (12 allocations: 288 bytes)
    # 2.78303613775333
end