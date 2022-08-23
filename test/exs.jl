using GaussianDispersion
using GaussianDispersion: potential_temp
using Plots

Rg = 8.31447
g = 9.80665
Mₐ = 0.028964
R = Rg / Mₐ

# 5.1 p.63
rawinsonde = [
    946 31
    938	26.8
    904	24.5
    900	24.1
    850	19.5
    800	14.6
    732	7.4
    702	5.5
    700	5.4
    666	4.1
    655	5.1
    633	3.8
    600	1
]

p_rw = rawinsonde[:,1] * 1e2
t_rw = rawinsonde[:,2] .+ 273.15

rho_rw = p_rw ./ (R .* t_rw)

zs = trapz(p_rw, rho_rw)

function trapz(p, rho, init = 0)
    fs = - 1 / g .* 0.5 .* (p[2:end] - p[1:end-1]) .* (1 ./ rho[1:end-1] .+ 1 ./ rho[2:end])
    ccat = vcat(init, init .+ fs...)
    cumsum(ccat)
end

tplot = plot(t_rw, zs, marker = :dot)

# 5.2 p.68
Γ = 0.00975 # Km-1
T0 = t_rw[1]
Tend = t_rw[1] - Γ * zs[end]

plot!(tplot, [T0, Tend], [zs[1], zs[end]], linestyle = :dash, color = :red)

# 5.3 p.69
Θ = potential_temp.(t_rw, p_rw, R, 1006)
plot(Θ, zs, marker = :dot)