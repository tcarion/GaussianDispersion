# GaussianDispersion

[![lifecycle](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tcarion.github.io/GaussianDispersion.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tcarion.github.io/GaussianDispersion.jl/dev/)
[![Build Status](https://github.com/tcarion/GaussianDispersion.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tcarion/GaussianDispersion.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/tcarion/GaussianDispersion.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/tcarion/GaussianDispersion.jl)

GaussianDispersion.jl is a Julia package for running Gaussian plume dispersion models. This package is under heavy development, it still lacks of a lot of features and is subjected to huge changes. If you want to contribute, don't hesitate to contact me or to pull request!

## Installation
The package is not yet on the official registry, and must be installed this way:
```julia
using Pkg; Pkg.add(url="https://github.com/tcarion/GaussianDispersion.jl")
```

# Quick start
It's currently possible to run a simple Gaussian plume model with the Pasquill and Gifford parametrization:

```julia
using GaussianDispersion

relpar = ReleaseParams(h = 12, Q = 100, u = 4)
plume = GaussianPlume(release = relpar)
plume.stabilities = Stabilities(:D)
plume.reflection = true

# Calculate the plume on a 3-D domain:
result = [plume(x,y,z) for x in range(0, 2000, 200), y in range(-300, 300, 100), z in range(0, 150, 150)]
```