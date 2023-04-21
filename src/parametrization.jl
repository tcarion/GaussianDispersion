abstract type AbstractGaussianModel end

# Models parametrization dispatch singletons
struct Aermod <: AbstractGaussianModel end
struct Calpuff  <: AbstractGaussianModel end

abstract type AbstractTerrain end

# Dispatch singletons for the type of terrains
struct Urban <: AbstractTerrain end
struct Rural <: AbstractTerrain end
struct Irrigated <: AbstractTerrain end
struct Water <: AbstractTerrain end
struct Snow <: AbstractTerrain end
struct Ice <: AbstractTerrain end
struct Forest <: AbstractTerrain end
struct Mountains <: AbstractTerrain end
struct Farmland <: AbstractTerrain end
struct Prairie <: AbstractTerrain end
struct Lawn <: AbstractTerrain end

Base.broadcastable(t::AbstractTerrain) = Ref(t)

# abstract type AbstractDispersionParameters end

# # struct DispersionCoefficientsTrait end
# struct BriggsDispersion{T, S} <: AbstractDispersionParameters 
#     terrain::T
#     stability::S
# end

# BriggsDispersion(terrain::AbstractTerrain, stab::PGStability) = BriggsDispersion(terrain, stab.class)
# # isdispersion(::Type{BriggsDispersion}) = true

# # dispersion_params(p::T) where T = dispersion_params(T, p)
# sigma_y(dispparams::AbstractDispersionParameters) = sigma_y(disp_coefs_functions(dispparams))
# sigma_z(dispparams::AbstractDispersionParameters) = sigma_z(disp_coefs_functions(dispparams))

"""
    AbstractDispersionFunctions
Define types for representing functions to compute the dispersion coefficients with respect to the downwind distance `x` in meter.
Those types must implement the `sigma_y` and `sigma_z` methods that both return a closure for calculating `σ_y` and `σ_z` respectively.
Instanciating those types give a functor that takes as argument the downwind distance and return the dispersion coefficients.
"""
abstract type AbstractDispersionFunctions end

_disp_fun_error(o) = error("You have not implemented the `sigma_y` or `sigma_z` functions for the object $(typeof(o))")
sigma_y(o::AbstractDispersionFunctions)::Function = _disp_fun_error(o)
sigma_z(o::AbstractDispersionFunctions)::Function = _disp_fun_error(o)

(funs::AbstractDispersionFunctions)(x::Number) = (sigma_y(funs)(x), sigma_z(funs)(x))