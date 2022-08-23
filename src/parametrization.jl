abstract type AbstractGaussianModels end

# Models parametrization dispatch singletons
struct Aermod <: AbstractGaussianModels end
struct Calpuff  <: AbstractGaussianModels end

abstract type AbstractTerrain end

# Dispatch singletons for the type of terrains
struct Urban <: AbstractTerrain end
struct Rural <: AbstractTerrain end
struct Irrigated <: AbstractTerrain end
struct Water <: AbstractTerrain end
struct Forest <: AbstractTerrain end

Base.broadcastable(t::AbstractTerrain) = Ref(t)