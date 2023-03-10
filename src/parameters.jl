"""
    $(TYPEDEF)
Structure related to the release conditions

    $(FIELDS)
"""
Base.@kwdef mutable struct ReleaseParams
    "Effective source height [m]"
    h::Real = 1.
    "Emission rate [g/s]"
    Q::Real = 1.
end

Base.@kwdef mutable struct DispersionCoefficients
    "Horizontal dispersion parameter [m]"
    y::Real = 2.15
    "Vertical dispersion parameter [m]"
    z::Real = 2.15
end

"""
    DispersionCoefficients(x::Real, terrain::AbstractTerrain, stability::StabilityClass)  

Return a `DispersionCoefficients` struct given the downwind direction `x`, the terrain and stability class
# Example
julia> DispersionCoefficients(100, Rural, A)
DispersionCoefficients(21.89081818461976, 20.0)
```
"""
# function DispersionCoefficients(x::Real, terrain::AbstractTerrain, stability::StabilityClass)
#     ycoef, zcoef = briggs_dispersion(typeof(terrain), stability)
#     sigma_y = disp_function(ycoef...)(x)
#     sigma_z = disp_function(zcoef...)(x)
#     DispersionCoefficients(sigma_y, sigma_z)
# end
# Base.:+(d1::DispersionCoefficients, d2::DispersionCoefficients) = DispersionCoefficients(d1.y + d2.y, d1.z + d2.z) 
# Base.:/(d1::DispersionCoefficients, x::Real) = DispersionCoefficients(d1.y / x, d1.z / x) 
# Base.:(==)(d1::DispersionCoefficients, d2::DispersionCoefficients) = d1.z == d2.z && d1.y == d2.y

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


struct BriggsFunctions <: AbstractDispersionFunctions
    σ_yx::Function
    σ_zx::Function
end

BriggsFunctions(a₁::T, b₁::T, c₁::T, a₂::T, b₂::T, c₂::T) where T = BriggsFunctions(_briggs_function(a₁, b₁, c₁), _briggs_function(a₂, b₂, c₂))

BriggsFunctions(t::AbstractTerrain, s::PGStability) = BriggsFunctions(t, s.class)
BriggsFunctions(t::AbstractTerrain, ss::Vector{PGStability}) = BriggsFunctions(t, [s.class for s in ss])

# This creates a new BriggsFunctions structure, modifying the encapsulated functions to be the mean of the dispersion coefficient values.
function Statistics.mean(funs::Vector{BriggsFunctions})
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

