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
    "parameters related to the atmospheric conditions"
    meteo::M = MeteoParams()
    "object representing the dispersion coefficients"
    dispcoefs::T = BriggsFunctions(Rural(), meteo.stability)
    "if ground reflection is considered"
    reflection::Bool = true

    _z_exp_terms::Vector{<:Real} = _construct_z_terms(release.h, meteo.hmix, reflection)
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
    concentration(y, z, wind, Q, coefs, plume._z_exp_terms)
end

# (plume::GaussianPlume)(x::Real, y::Real, z::Real) = concentration(x, y, z, plume.release, plume.dispcoefs; reflection = plume.reflection)

mixing_layer_additional_term(hmix, h) = vcat([[-h - 2*(i+1)*hmix, h + 2*(i+1)*hmix, -h + 2*(i+1)*hmix, h - 2*(i+1)*hmix] for i in 0:9]...)
downwind_mass(Q, u) = Q / u
gaussian_exp(x, p1, p2) = exp(-0.5 * ((x + p1) / p2)^2)
crosswind_conc(y, sigma_y) = 1 / sqrt(2pi) / sigma_y * gaussian_exp(y, 0, sigma_y)
vertical_conc(z, sigma_z, h) = 1 / sqrt(2pi) / sigma_z * sum(gaussian_exp.(z, h, sigma_z))
function concentration(y, z, u, Q, disp::DispersionCoefficients, z_exp_terms)
    sigy = disp.y
    sigz = disp.z
    downwind_mass(Q, u) * crosswind_conc(y, sigy) * vertical_conc(z, sigz, z_exp_terms)
end

function _construct_z_terms(h, hmix, reflection) 
    z_exp_terms = [-h]
    z_exp_terms = reflection ? [h, -h] : z_exp_terms
    z_exp_terms = !isnothing(hmix) ? vcat(z_exp_terms, mixing_layer_additional_term(hmix, h)) : z_exp_terms
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
