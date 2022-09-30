
@enum StabilityClass A B C D E F

const StabilityClasses = AbstractSet{StabilityClass}

"""
    $(TYPEDSIGNATURES)
Define the stability classes for the model according to Pasquill and Gifford (Pasquill, 1961; Gifford, 1961).
If multiple classes are defined, the average result for each class is considered.

# Example
```jl
julia> Stabilities(:A, :B)
Set{StabilityClass} with 2 elements:
  GaussianDispersion.A
  GaussianDispersion.B
```
"""
function Stabilities(stabs::Vararg{Union{String, Symbol}}) :: StabilityClasses
    Set(tostab(stabs...))
end

function tostab(stabs::Vararg{Symbol})
    [getproperty(@__MODULE__, s) for s in collect(stabs)]
end

function tostab(stabs::Vararg{AbstractString})
    sstabs = Symbol.(stabs)
    tostab(sstabs...)
end

"""
    briggs_dispersion(terrain::AbstractTerrain, stability::StabilityClass)

Given the `terrain` and the `stability`, return the three coefficients a, b, c for σ_y and σ_z respectively.

# Example
```jl
julia> briggs_dispersion(Rural, A)
((a = 0.22, b = 0.0001, c = -0.5), (a = 0.2, b = 0, c = 1.0))
```
"""
function briggs_dispersion(::Type{Rural}, stability::StabilityClass)
    if stability == A
        sigycol = 1
        sigzcol = 7
    elseif stability == B
        sigycol = 2
        sigzcol = 8
    elseif stability == C
        sigycol = 3
        sigzcol = 9
    elseif stability == D 
        sigycol = 4
        sigzcol = 10
    elseif stability == E
        sigycol = 5
        sigzcol = 11
    elseif stability == F
        sigycol = 6
        sigzcol = 12
    end
    briggs_dispersion(sigycol, sigzcol)
end

function briggs_dispersion(::Type{Urban}, stability::StabilityClass)
    if stability == A || stability == B
        sigycol = 13
        sigzcol = 17
    elseif stability == C
        sigycol = 14
        sigzcol = 18
    elseif stability == D
        sigycol = 15
        sigzcol = 19
    elseif stability == E || stability == F
        sigycol = 16
        sigzcol = 20
    end
    briggs_dispersion(sigycol, sigzcol)
end

briggs_dispersion(sigycol::Int, sigzcol::Int) = (Tuple(BRIGGS_COEFS[:, sigycol]), Tuple(BRIGGS_COEFS[:, sigzcol]))
"""
    disp_function(a, b, c) 

Return the function for calculating the dispersion parameters with respect to x given the coefficients of the equation ax(1 + bx)^c
```
"""
disp_function(a, b, c) = x -> a * x * (1 + b*x)^c

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
    "Wind speed [m/s]"
    u::Real = 5.
end

Base.@kwdef mutable struct DispersionParams
    "Horizontal dispersion parameter [m]"
    y::Real = 2.15
    "Vertical dispersion parameter [m]"
    z::Real = 2.15
end

"""
    DispersionParams(x::Real, terrain::AbstractTerrain, stability::StabilityClass)  

Return a `DispersionParams` struct given the downwind direction `x`, the terrain and stability class
# Example
julia> DispersionParams(100, Rural, A)
DispersionParams(21.89081818461976, 20.0)
```
"""
function DispersionParams(x::Real, terrain::AbstractTerrain, stability::StabilityClass)
    ycoef, zcoef = briggs_dispersion(typeof(terrain), stability)
    sigma_y = disp_function(ycoef...)(x)
    sigma_z = disp_function(zcoef...)(x)
    DispersionParams(sigma_y, sigma_z)
end
Base.:+(d1::DispersionParams, d2::DispersionParams) = DispersionParams(d1.y + d2.y, d1.z + d2.z) 
Base.:/(d1::DispersionParams, x::Real) = DispersionParams(d1.y / x, d1.z / x) 
Base.:(==)(d1::DispersionParams, d2::DispersionParams) = d1.z == d2.z && d1.y == d2.y

"""
    $(TYPEDEF)
Values for the heat balance parametrization.

- 'albedo' : 
- `C_g`: parameter for the estimation of the heat flux to the ground in the equation ``q\\_G = C\\_G R\\_n``. 
Typical values are `0.05-0.25` for rural, `0.25-0.3` for urban and `0.1` for grasslands (Scire et al., 2000)
"""
Base.@kwdef mutable struct HeatBalanceParams
    C_g::Real = 0.1
    albedo::Real = 0.18
    B::Real = 1.
end

HeatBalanceParams(::Type{Aermod}) = HeatBalanceParams(
    C_g = 0.1
)

HeatBalanceParams(::Type{Calpuff}, terrain::Type{<:AbstractTerrain}) = HeatBalanceParams(;
    CALPUFF_HEAT_PARAMS[Symbol(terrain)]...
)

# p.78
CALPUFF_HEAT_PARAMS = (
    Urban       = (albedo = 0.18, B = 1.5),
    Rural       = (albedo = 0.20, B = 1.0),
    Irrigated   = (albedo = 0.15, B = 0.5),
    Water       = (albedo = 0.10, B = 0.0),
    Forest      = (albedo = 0.10, B = 1.0),
)

function sensible_heat_flux(hbp::HeatBalanceParams, cloud_cover, T_surf, sol_elev)
    @unpack C_g, albedo, B = hbp
    R = solar_energy_flux(sol_elev, cloud_cover)
    net_rad = net_radiation(R, albedo, cloud_cover, T_surf)
    sensible_heat_flux(C_g, net_rad, B)
end