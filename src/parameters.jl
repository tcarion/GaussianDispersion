"""
    $(TYPEDEF)
Structure related to the release conditions

    $(FIELDS)
"""
Base.@kwdef struct ReleaseParams{T}
    "Effective source height [m]"
    h::Real = 1.
    "Emission rate [g/s]"
    Q::Real = 1.
    "Gas temperature out of the stack [K]"
    Ts::T = nothing
end

Base.@kwdef struct DispersionCoefficients
    "Horizontal dispersion parameter [m]"
    y::Real = 2.15
    "Vertical dispersion parameter [m]"
    z::Real = 2.15
end

"""
    $(SIGNATURES)
Buoyancy flux parameter according to Briggs (1968).
"""
briggs_buoyancy_flux(gas_density, air_density, stack_radius, gas_velocity) = (1 - gas_density / air_density) * GRAVITY * stack_radius^2 * gas_velocity

"""
    $(SIGNATURES)
Plume rise `Δh` [m] according to Briggs (1975) parametrization. `x` is the downwind distance in meter and `u` the downwind velocity in meter/second.
"""
function plume_rise(flux_param, x, u) 
    if flux_param <= 55.
        xf = 49. * flux_param^(5/8)
    else
        xf = 119. * flux_param^(2/5)
    end
    1.6 * flux_param^(1/3) * min(x, xf)^(2/3) / u
end

"""
    $(TYPEDEF)
Values for the heat balance parametrization.

- 'albedo' : 
- `C_g`: parameter for the estimation of the heat flux to the ground in the equation ``q\\_G = C\\_G R\\_n``. 
Typical values are `0.05-0.25` for rural, `0.25-0.3` for urban and `0.1` for grasslands (Scire et al., 2000)
""" 
Base.@kwdef struct HeatBalanceParams{C, A, B}
    C_g::C = Param(0.1, bounds = (0.05, 0.3))
    albedo::A = Param(0.18, bounds = (0.1, 0.18))
    "Bowen ratio"
    Bowen::B = Param(1., bounds = (0., 10.))
end

HeatBalanceParams(::Type{Aermod}) = HeatBalanceParams(
    C_g = 0.1
)

HeatBalanceParams(::Type{Calpuff}, terrain::Type{<:AbstractTerrain}) = HeatBalanceParams(;
    CALPUFF_HEAT_PARAMS[Symbol(terrain)]...
)

# p.78
CALPUFF_HEAT_PARAMS = (
    Urban       = (albedo = 0.18, Bowen = 1.5),
    Rural       = (albedo = 0.20, Bowen = 1.0),
    Irrigated   = (albedo = 0.15, Bowen = 0.5),
    Water       = (albedo = 0.10, Bowen = 0.0),
    Forest      = (albedo = 0.10, Bowen = 1.0),
)

function sensible_heat_flux(hbp::HeatBalanceParams, cloud_cover, T_surf, sol_elev)
    @unpack C_g, albedo, B = hbp
    R = solar_energy_flux(sol_elev, cloud_cover)
    net_rad = net_radiation(R, albedo, cloud_cover, T_surf)
    sensible_heat_flux(C_g, net_rad, B)
end


struct BriggsFunctions{F1<:Function, F2<:Function} <: AbstractDispersionFunctions
    σ_yx::F1
    σ_zx::F2
end

BriggsFunctions(a₁::T, b₁::T, c₁::T, a₂::T, b₂::T, c₂::T) where T = BriggsFunctions(_briggs_function(a₁, b₁, c₁), _briggs_function(a₂, b₂, c₂))

BriggsFunctions(t::AbstractTerrain, s::PGStability) = BriggsFunctions(t, s.class)
BriggsFunctions(t::AbstractTerrain, ss::Vector{PGStability}) = BriggsFunctions(t, [s.class for s in ss])

# This creates a new BriggsFunctions structure, modifying the encapsulated functions to be the mean of the dispersion coefficient values.
function Statistics.mean(funs::Vector{<:BriggsFunctions})
    BriggsFunctions(
        x -> mean(x .|> sigma_y.(funs)),
        x -> mean(x .|> sigma_z.(funs)),
    )
end

# If several stability classes are given, we consider the mean of the coefficients.
BriggsFunctions(t::AbstractTerrain, s::Vector{PGStabilityClass}) = mean(BriggsFunctions.(t, s))

BriggsFunctions(::Rural, ::PGVeryUnstable)      = BriggsFunctions(BRIGGS_COEFS[:, 1]..., BRIGGS_COEFS[:, 7]...)
BriggsFunctions(::Rural, ::PGUnstable)          = BriggsFunctions(BRIGGS_COEFS[:, 2]..., BRIGGS_COEFS[:, 8]...)
BriggsFunctions(::Rural, ::PGSlightlyUnstable)  = BriggsFunctions(BRIGGS_COEFS[:, 3]..., BRIGGS_COEFS[:, 9]...)
BriggsFunctions(::Rural, ::PGNeutral)           = BriggsFunctions(BRIGGS_COEFS[:, 4]..., BRIGGS_COEFS[:, 10]...)
BriggsFunctions(::Rural, ::PGSlightlyStable)    = BriggsFunctions(BRIGGS_COEFS[:, 5]..., BRIGGS_COEFS[:, 11]...)
BriggsFunctions(::Rural, ::PGStable)            = BriggsFunctions(BRIGGS_COEFS[:, 6]..., BRIGGS_COEFS[:, 12]...)
BriggsFunctions(::Urban, ::Union{PGVeryUnstable, PGUnstable}) = BriggsFunctions(BRIGGS_COEFS[:, 13]..., BRIGGS_COEFS[:, 17]...)
BriggsFunctions(::Urban, ::PGSlightlyUnstable)                = BriggsFunctions(BRIGGS_COEFS[:, 14]..., BRIGGS_COEFS[:, 18]...)
BriggsFunctions(::Urban, ::PGNeutral)                         = BriggsFunctions(BRIGGS_COEFS[:, 15]..., BRIGGS_COEFS[:, 19]...)
BriggsFunctions(::Urban, ::Union{PGSlightlyStable, PGStable}) = BriggsFunctions(BRIGGS_COEFS[:, 16]..., BRIGGS_COEFS[:, 20]...)


sigma_y(relation::BriggsFunctions) = relation.σ_yx
sigma_z(relation::BriggsFunctions) = relation.σ_zx


_briggs_function(a::T, b::T, c::T) where T = x -> a * x * (1 + b*x)^c

DispersionCoefficients(briggs::BriggsFunctions, x) = DispersionCoefficients(briggs(x)...) 

