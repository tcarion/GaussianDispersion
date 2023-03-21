using GaussianDispersion

const GD = GaussianDispersion

const coefs = [
#   A      α      B      β
    0.235  0.796  0.311  0.711
    0.297  0.796  0.382  0.711
    0.418  0.796  0.52   0.711
    0.586  0.796  0.7    0.711
    0.826  0.796  0.95   0.711
    0.946  0.796  1.321  0.711
    1.043  0.698  0.819  0.669
]

abstract type BMStabilityClass end

struct BMVeryStable <: BMStabilityClass end
struct BMStable <: BMStabilityClass end
struct BMNeutral <: BMStabilityClass end
struct BMSlightlyUnstable <: BMStabilityClass end
struct BMUnstable <: BMStabilityClass end
struct BMVeryUnstable <: BMStabilityClass end
struct BMHighWind <: BMStabilityClass end

# We define a struct implementing the `AbstractDispersionFunctions` class.
struct BMFunctions{F1<:Function, F2<:Function} <: GD.AbstractDispersionFunctions
    f1::F1
    f2::F2
end

# We define a constructor that builds the dispersion functions given the 4 coefficients from Bultinck.
BMFunctions(A::Number, α::Number, B::Number, β::Number) = BMFunctions(_bm_function(A, α), _bm_function(B, β)) 
_bm_function(A, α) = x -> A * x^α

# Multiple functions are built given the stability class.
BMFunctions(::BMVeryStable) =       BMFunctions(coefs[1, :]...)
BMFunctions(::BMStable) =           BMFunctions(coefs[2, :]...)
BMFunctions(::BMNeutral) =          BMFunctions(coefs[3, :]...)
BMFunctions(::BMSlightlyUnstable) = BMFunctions(coefs[4, :]...)
BMFunctions(::BMUnstable) =         BMFunctions(coefs[5, :]...)
BMFunctions(::BMVeryUnstable) =     BMFunctions(coefs[6, :]...)
BMFunctions(::BMHighWind) =         BMFunctions(coefs[7, :]...)

# We extend both methods providing the dispersion functions.
GD.sigma_y(bm::BMFunctions) = bm.f1
GD.sigma_z(bm::BMFunctions) = bm.f2

bm = BMFunctions(BMStable())

bm(5.)

relpar = ReleaseParams(h = 75, Q = 100)

plume_briggs_veryunstable_rural = GaussianPlume(release = relpar, meteo = MeteoParams(stability = PGStability(:A)))
plume_briggs_veryunstable = GaussianPlume(release = relpar, dispcoefs = GD.BriggsFunctions(Urban(), PGStability(:A)))
plume_briggs_stable = GaussianPlume(release = relpar, dispcoefs = GD.BriggsFunctions(Urban(), PGStability(:E)))
plume_briggs_stable_rural = GaussianPlume(release = relpar, dispcoefs = GD.BriggsFunctions(Rural(), PGStability(:E)))

plume_verystable = GaussianPlume(release = relpar, dispcoefs = BMFunctions(BMVeryStable()))
plume_stable = GaussianPlume(release = relpar, dispcoefs = BMFunctions(BMStable()))
plume_veryunstable = GaussianPlume(release = relpar, dispcoefs = BMFunctions(BMVeryUnstable()))

x, y, z = 150., 0., 0.
plume_briggs_stable(x, y, z)
plume_stable(x, y, z)
plume_veryunstable(x, y, z)


##
xs = range(0., 1000., length = 200)
plot(xs, plume_briggs_veryunstable_rural.(xs, y, z), label ="briggs very unstable", linestyle=:dash)
plot!(xs, plume_briggs_unstable_urban.(xs, y, z), label ="briggs very unstable rural", linestyle=:dash)
plot!(xs, plume_unstable.(xs, y, z), label = "bultinck very unstable", linestyle=:dash)
# plot!(xs, plume_verystable.(xs, y, z), label = "bultinck verystable")
plot!(xs, plume_briggs_stable.(xs, y, z), label ="briggs stable")
plot!(xs, plume_briggs_stable_rural.(xs, y, z), label ="briggs stable rural")
plot!(xs, plume_stable.(xs, y, z), label = "bultinck stable")
##