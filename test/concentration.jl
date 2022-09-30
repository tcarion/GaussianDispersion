using GaussianDispersion

# Example from reference book "Air Dispersion Modeling - Foundations and Applications, Alex De Visscher, p.24-25"
@testset "concentration" begin
    # Example 2.1
    u = 7.
    h = 75 + 15
    stabs = pasquill_gifford(Moderate, u)
    relpar = ReleaseParams(h = h, Q = 100, u = u)
    plume = GaussianPlume(release = relpar)
    plume.stabilities = stabs
    cground1 = plume(1500, 0, 0) * 1e6
    ccross = plume(1500, 100, 0) * 1e6

    @test cground1 ≈ 160.3 atol=1e-1
    @test ccross ≈ 107.5 atol=1e-1

    h = 30
    relpar = ReleaseParams(h = h, Q = 5, u = 2)

    plume = GaussianPlume(release = relpar)
    
    plume.stabilities = Stabilities(:D)
    plume.reflection = true

    cground2 = plume(1000, 0, 0)
    @test cground2 ≈ 201.139895390322 * 1e-6
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
    plume = GaussianPlume()
    plume.stabilities = Stabilities(:A, :B)
    cboth= plume(1000, 0, 0)

    plume.stabilities = Stabilities(:A)
    ca= plume(1000, 0, 0)

    plume.stabilities = Stabilities(:B)
    cb= plume(1000, 0, 0)
end