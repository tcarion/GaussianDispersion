
abstract type AbstractMeteoParams end

abstract type AbstractSkyCondition end

abstract type AbstractStability end

"""
    PGStabilityClass
Abstract type of Pasquill-Gifford stability classes. The concrete types start with `PG` followed by the letter from Pasquill-Gifford classes (A, B, C, D, E, F).

# Reference
(Pasquill, 1961; Gifford, 1961).
"""
abstract type PGStabilityClass end

struct PGVeryUnstable <: PGStabilityClass end
struct PGUnstable <: PGStabilityClass end
struct PGSlightlyUnstable <: PGStabilityClass end
struct PGNeutral <: PGStabilityClass end
struct PGSlightlyStable <: PGStabilityClass end
struct PGStable <: PGStabilityClass end

const MAP_PG_CLASSES = (
    A = PGVeryUnstable(),
    B = PGUnstable(),
    C = PGSlightlyUnstable(),
    D = PGNeutral(),
    E = PGSlightlyStable(),
    F = PGStable(),
)

"""
    PGStability <: AbstractPGStability
Define the stability classes for the model according to Pasquill and Gifford.
#! If multiple classes are defined, the average result for each class is considered.

# Example
```jl-repl
julia> pgstab = PGStability(:A)
PGStability(GaussianDispersion.PGVeryUnstable())
```
"""
struct PGStability <: AbstractStability
    class::PGStabilityClass
end

PGStability(letter::Symbol) = PGStability(MAP_PG_CLASSES[letter])

# function PGStability(classes::Vararg{Symbol})
#     pgclasses = map(_symbol_to_pg, classes)
#     PGStability(Set(pgclasses))
# end

# function _symbol_to_pg(s::Symbol)
#     object_name = Meta.parse("PG$s")
#     eval(:($object_name()))
# end

struct Strong <: AbstractSkyCondition end
struct Moderate <: AbstractSkyCondition end
struct Slight <: AbstractSkyCondition end
struct Cloudy <: AbstractSkyCondition end
struct Clear <: AbstractSkyCondition end

"""
    pasquill_gifford(criteria::AbstractSkyCondition, windspeed::Real)
Return a list of stability classes according to the Pasquill Gifford criteria

# Example
julia> pasquill_gifford(Moderate, 5.5)
Set{StabilityClass} with 2 elements:
 C
 D
"""
function pasquill_gifford(criteria::AbstractSkyCondition, windspeed::Real)
    if windspeed < 2
        if criteria == Strong()
            Set([PGVeryUnstable()])
        elseif criteria == Moderate()
            Set([PGVeryUnstable(), PGUnstable()])
        elseif criteria == Slight()
            Set([PGUnstable()])
        elseif criteria == Cloudy()
            Set([PGSlightlyStable()])
        elseif criteria == Clear()
            Set([PGStable()])
        end
    elseif 2 <= windspeed < 3
        if criteria == Strong()
            Set([PGVeryUnstable(), PGUnstable()])
        elseif criteria == Moderate()
            Set([PGUnstable()])
        elseif criteria == Slight()
            Set([PGSlightlyUnstable()])
        elseif criteria == Cloudy()
            Set([PGSlightlyStable()])
        elseif criteria == Clear()
            Set([PGStable()])
        end
    elseif 3 <= windspeed < 5
        if criteria == Strong()
            Set([PGUnstable()])
        elseif criteria == Moderate()
            Set([PGUnstable(), PGSlightlyUnstable()])
        elseif criteria == Slight()
            Set([PGSlightlyUnstable()])
        elseif criteria == Cloudy()
            Set([PGNeutral()])
        elseif criteria == Clear()
            Set([PGSlightlyStable()])
        end
    elseif 5 <= windspeed <= 6
        if criteria == Strong()
            Set([PGSlightlyUnstable()])
        elseif criteria == Moderate()
            Set([PGSlightlyUnstable(), PGNeutral()])
        elseif criteria == Slight()
            Set([PGNeutral()])
        elseif criteria == Cloudy()
            Set([PGNeutral()])
        elseif criteria == Clear()
            Set([PGNeutral()])
        end
    else
        if criteria == Strong()
            Set([PGSlightlyUnstable()])
        elseif criteria == Moderate()
            Set([PGNeutral()])
        elseif criteria == Slight()
            Set([PGNeutral()])
        elseif criteria == Cloudy()
            Set([PGNeutral()])
        elseif criteria == Clear()
            Set([PGNeutral()])
        end
    end
end
Base.@kwdef struct MeteoParams{W, S, M} <: AbstractMeteoParams
    wind::W = 5.
    stability::S = PGStability(:A)
    "height of the mixing layer [m]"
    hmix::M = nothing
end

"""
    $(TYPEDSIGNATURES)

Calculate Obukhov scale height from surface meteorological data and sensible heat flux.
# Arguments                                           
- ps:      surface pressure [Pa]                  
- ts:      surface temperature [K]                
- td:      surface dew point [K]                  
- t:       temperature first model level [K]
- p :      pressure first model level [Pa]
- u_star: scale velocity [m/s]                   
- q:       surface sensible heat flux [W/m2]      
"""
function obukhov(ps, ts, td, t, p, u_star, q, c_p = 1004.5, R_gas = 287.05, karman = 0.4, g = 9.81)
    e = saturation_pressure(td)
    tv = virtual_temp(ps, t, e)
    rho = p / (R_gas * tv)
    θ = potential_temp(t, p, c_p, R_gas)
    θ_star = friction_temp(q, rho, u_star, c_p)
    u_star^2 * ts / (karman * g * θ_star)
end
"""
    $(TYPEDSIGNATURES)

Calculate the saturation vapor pressure from Arden Buck equations (Buck, 1996). `t` is in K, result in Pa.
"""
function saturation_pressure(t)
    t = t - 273.15
    if t >= 0.
        r = 6.1121 * exp( (18.678 - t / 234.5) * t / (257.14 + t) )
    else
        r = 6.1115 * exp( (23.036 - t / 333.7) * t / (279.82 + t) ) 
    end
    r * 100.
end

"""
    $(TYPEDSIGNATURES)

Calculate the virtual temperature from the surface pressure `ps` [Pa], the temperature `t` [K], and the vapor pressure `e` [Pa].
"""
virtual_temp(ps, t, e) = t / (1 - 0.378 * e / ps )

"""
    $(TYPEDSIGNATURES)

Calculate the virtual temperature from the temperature `t` [K] and the specific humidity `w` [kg/kg].
"""
virtual_temp(t, w) = t * (1 + 0.608 * w )

"""
    $(TYPEDSIGNATURES)

Calculate the potential temperature.

# Arguments                                           
- t:       temperature [K]                
- p:       pressure [Pa]                
- R_gas:   individual gas constant for dry air [J/kg/K]
- c_p:     specific heat capacities of dry air [J/kg/K]                
"""
potential_temp(t, p, akap) = t * (1e5 / p)^(akap)
potential_temp(t, p, R_gas, c_p) = potential_temp(t, p, R_gas / c_p)


"""
    $(TYPEDSIGNATURES)

Calculate the friction temperature from the surface pressure `p` [Pa], temperature `t` [K], dewpoint temp `td` [K] and stress `stress` [N/m²].

# Arguments                                              
- q:       surface sensible heat flux [W/m²]
- rho:     air density [kg/m³]
- u_star:   scale velocity [m/s]
- c_p:     specific heat capacities of dry air [J/kg/K]
"""
friction_temp(q, rho, u_star, c_p) = - q / (rho * c_p * u_star)

"""
    $(TYPEDSIGNATURES)

Calculate the friction velocity from the surface pressure `p` [Pa], temperature `t` [K], dewpoint temp `td` [K] and stress `stress` [N/m²].
"""
function friction_velocity(p, t, td, stress, R_gas)
    # vapor pressure
    e = saturation_pressure(td)
    tv = virtual_temp(p, t, e)
    rho_moist = p / (R_gas * tv)
    sqrt(abs(stress) / rho_moist)
end

"""
    $(TYPEDSIGNATURES)

Get height from pressure level in a hydrostatic atmosphere.
"""
height_hs(p, ps, R_gas, T, g) = -R_gas * T / g * log(p/ps)

"""
    $(TYPEDSIGNATURES)

Get height from pressure level in the international standard atmosphere.
"""
height_isa(p, R_gas, g) = height_hs(p, 101325, R_gas, 288, g)

## Empirical parametrizations

### Related to heat balance (≈ p.78)
"""
    $(TYPEDSIGNATURES)

Solar radiation energy flux W m^-2 given the solar elevation `Φ` and `n` the fractional cloud cover.
Ref: Kasten and Czeplak, 1980; Holtslag and van Ulden, 1983.
"""
solar_energy_flux(Φ, n) = (990 * sind(Φ) - 30) * (1 - 0.75 * n^3.4)

"""
    $(TYPEDSIGNATURES)

Net radiation energy flux received by the earth surface [W m^-2] with `R` the solar flux ([`solar_energy_flux`](@ref))
"""
net_radiation(R, albedo, cloud_cover, T_surf) = ( R*(1-albedo) + 60 * cloud_cover - 5.67 * 1e-8 * T_surf^4 + 5.31e-13 * T_surf^6 ) / 1.12

"""
    $(TYPEDSIGNATURES)

Sensible heat flux to the surface [W m^-2]. `net_rad` can be calculated with [`net_radiation`](@ref)
"""
sensible_heat_flux(C_g, net_rad, B) = (1 - C_g) * net_rad / (1 + 1 / B)

"""
    $(TYPEDSIGNATURES)

Latitude [°] where the sun is in the zenith at the solar noon. `day` is the day of the year
"""
solar_declination(day) = 23.45 * cosd(360 * (day + 10) / 365.25)

# The equation of the reference book (p.558) gives bad results
# solar_declination(day) = 23.45 * sind(360 * (day - 284) / 365.25)

"""
    $(TYPEDSIGNATURES)

`t` the time in hour and `t0` the time of the solar noon.
"""
hour_angle(t, t0 = 12) = 15 * (t - t0)

"""
    $(TYPEDSIGNATURES)

Solar elevation [°] given the latitude, the solar declination and the hour angle.
"""
solar_elevation(lat, sol_dec, h_angle) = asind(sind(lat) * sind(sol_dec) + cosd(lat) * cosd(sol_dec) * cosd(h_angle))

"""
    $(TYPEDSIGNATURES)

Albedo given `albedo_90` the albedo at 90° and `ϕ` the solar elevation
"""
albedo_elevation(albedo_90, ϕ) = albedo_90 + (1 - albedo_90) * exp(-0.1 * ϕ - 0.5 * (1 - albedo_90)^2)