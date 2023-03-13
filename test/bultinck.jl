using GaussianDispersion

const GD = GaussianDispersion

coefs = [
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

struct BMFunctions{F1<:Function, F2<:Function} <: GD.AbstractDispersionFunctions
    f1::F1
    f2::F2
end

BMFunctions(::BMVeryStable) =       BMFunctions(coefs[1, :]...)
BMFunctions(::BMStable) =           BMFunctions(coefs[2, :]...)
BMFunctions(::BMNeutral) =          BMFunctions(coefs[3, :]...)
BMFunctions(::BMSlightlyUnstable) = BMFunctions(coefs[4, :]...)
BMFunctions(::BMUnstable) =         BMFunctions(coefs[5, :]...)
BMFunctions(::BMVeryUnstable) =     BMFunctions(coefs[6, :]...)
BMFunctions(::BMHighWind) =         BMFunctions(coefs[7, :]...)

BMFunctions(A::Number, α::Number, B::Number, β::Number) = BMFunctions(_bm_function(A, α), _bm_function(B, β)) 

GD.sigma_y(bm::BMFunctions) = bm.f1
GD.sigma_z(bm::BMFunctions) = bm.f2

_bm_function(A, α) = x -> A * x^α

bm = BMFunctions(BMStable())

bm(5.)