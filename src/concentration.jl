"""
    AbstractSkyCondition
Critera for determining the stability class according to sky conditions (Pasquill, 1961; Gifford, 1961)
`Strong`, `Moderate` and `Slight` are for incoming solar radiation (for daytime)
`Cloudy` and `Clear` are for cloudiness (for nighttime)
"""
# @enum AbstractSkyCondition Strong Moderate Slight Cloudy Clear

const DispersionParamsTypes = Union{DispersionCoefficients, AbstractDispersionFunctions}
"""
    $(TYPEDEF)
Structure related to the release conditions

    $(FIELDS)
"""
Base.@kwdef struct GaussianPlume{M, T}
    "release parameters"
    release::ReleaseParams = ReleaseParams()
    "terrain (Urban/Rural)"
    terrain::AbstractTerrain = Rural()
    "parameters related to the atmospheric conditions"
    meteo::M = MeteoParams()
    "object representing the dispersion coefficients"
    dispcoefs::T = BriggsFunctions(terrain, meteo.stability)
    "if ground reflection is considered"
    reflection::Bool = true
end

"""
    (plume::GaussianPlume)(x::Real, y::Real, z::Real)
Return the concentration in [g m^{-3}] at some point given the conditions stated in `plume`.
- `x` is the downwind distance from the point source.
- `y` is the horizontal distance perpendicular to the wind.
- `z` is the vertical distance from the point source.
"""
function (plume::GaussianPlume{M, <:AbstractDispersionFunctions})(x::Real, y::Real, z::Real) where M
    @unpack h, Q = plume.release
    @unpack wind = plume.meteo

    coefs = DispersionCoefficients(plume.dispcoefs(x)...)
    concentration(y, z, wind, h, Q, coefs; reflection = plume.reflection)
end

# (plume::GaussianPlume)(x::Real, y::Real, z::Real) = concentration(x, y, z, plume.release, plume.dispcoefs; reflection = plume.reflection)

"""
    $(SIGNATURES)
Buoyancy flux parameter according to Briggs (1968).
"""
buoyancy_flux(gas_density, air_density, stack_radius, gas_velocity) = (1 - gas_density / air_density) * GRAVITY * stack_radius^2 * gas_velocity

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

downwind_mass(Q, u) = Q / u
gaussian_exp(x, p1, p2) = exp(-0.5 * ((x + p1) / p2)^2)
crosswind_conc(y, sigma_y) = 1 / sqrt(2pi) / sigma_y * gaussian_exp(y, 0, sigma_y)
vertical_conc(z, sigma_z, h) = 1 / sqrt(2pi) / sigma_z * sum(gaussian_exp.(z, h, sigma_z))
function concentration(y, z, u, h, Q, disp::DispersionCoefficients; reflection = true)
    sigma_y = disp.y
    sigma_z = disp.z
    hs = reflection ? [h, -h] : -h
    downwind_mass(Q, u) * crosswind_conc(y, sigma_y) * vertical_conc(z, sigma_z, hs)
end

# function concentration(x::Real, y::Real, z::Real, release::ReleaseParams, coefs_fun::AbstractDispersionFunctions; reflection = true)
#     disp = DispersionCoefficients(sigma_y(coefs_fun)(x), sigma_z(coefs_fun)(x))
#     # If more than 1 stability, we take the average of the dispersion parameters for each class
#     # disp = sum(disps) / length(disps)
#     concentration(y, z, disp, release; reflection)
# end

# """
#     $(TYPEDSIGNATURES)
# Return the concentration in [g m^{-3}] at some point given the conditions stated in `params`.
# - `x` is the downwind distance from the point source.
# - `y` is the horizontal distance perpendicular to the wind.
# - `z` is the vertical distance from the point source.
# """
# concentration(x::Real, y::Real, z::Real, params::GaussianPlume) = 
#     concentration(x, y, z, params.release, params.terrain, params.stabilities; reflection = params.reflection)
