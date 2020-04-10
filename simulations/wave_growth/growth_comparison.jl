include("visualization.jl")
include("analysis.jl")
include("files.jl")

using Oceananigans, Printf

import Oceananigans: g_Earth

png_filepath = "growth_intro_plot.png"

wave_filename = "data/growth_wave_forced_convecting_spinup_ramp_a1.5_Nh256_Nz192_averages.jld2"
wind_filename = "data/growth_wind_forced_convecting_spinup_ramp_a1.5_Nh256_Nz192_averages.jld2"
both_filename = "data/growth_wind_with_waves_convecting_spinup_ramp_a1.5_Nh256_Nz192_averages.jld2"

aˢʷ, kˢʷ = get_surface_wave_parameters(wave_filename)
T = 4hour
Qᵘ₀ = get_parameter(wave_filename, "boundary_conditions", "Qᵘ₀")
  f = get_parameter(wave_filename, "coriolis", "f")

if Qᵘ₀ === nothing
    Qᵘ₀ = -5e-6
end

 Uˢ = aˢʷ^2 * kˢʷ * sqrt(g_Earth * kˢʷ)
Qᵘ₁ = -aˢʷ^2 * sqrt(g_Earth * kˢʷ) / 2T
  τ = abs(Qᵘ₁)

  u★ = sqrt(τ)

grid = get_grid(wave_filename)

t_w, U_w, V_w, B_w, Bz_w, u²_w, v²_w, w²_w = get_averages(wave_filename)
h_w, ∫U_w, ∫V_w, ∫u²_w, ∫v²_w, ∫w²_w = timeseries_from_averages(grid, U_w, V_w, B_w, u²_w, v²_w, w²_w)

t_τ, U_τ, V_τ, B_τ, Bz_τ, u²_τ, v²_τ, w²_τ = get_averages(wind_filename)
h_τ, ∫U_τ, ∫V_τ, ∫u²_τ, ∫v²_τ, ∫w²_τ = timeseries_from_averages(grid, U_τ, V_τ, B_τ, u²_τ, v²_τ, w²_τ)

t_b, U_b, V_b, B_b, Bz_b, u²_b, v²_b, w²_b = get_averages(both_filename)
h_b, ∫U_b, ∫V_b, ∫u²_b, ∫v²_b, ∫w²_b = timeseries_from_averages(grid, U_b, V_b, B_b, u²_b, v²_b, w²_b)

ramp(t, δ) = 1 - exp(-t^2 / (2δ^2))
∂t_ramp(t, δ) = exp(-t^2 / (2δ^2)) * t / δ^2

Nt = length(t_w)
Qᵘ(t) = Qᵘ₁ * T * ∂t_ramp(t, T) + Qᵘ₀
τmax = maximum(abs.(Qᵘ.(t_w)))

close("all")

#####
##### Final profiles
#####

fig, axs = subplots(ncols=5, sharey=true, figsize=(10, 4))

lw₁ = 2
lw₂ = 1.5
α₁ = 0.6
α₂ = 1.0

for (j, i) in enumerate((1, 201))

    if j == 1
        w_lbl = "Initial condition"
        τ_lbl = ""
        b_lbl = ""
    else
        w_lbl = @sprintf("Wave forced, \$ t = %d \\pi/f \$", f*t_w[i]/π)
        τ_lbl = @sprintf("Wind forced, \$ t = %d \\pi/f\$", f*t_w[i]/π)
        b_lbl = @sprintf("Wind forced beneath steady waves, \$ t = %d \\pi/f \$", f*t_w[i]/π)
    end

    sca(axs[1])
    plot(B_w[i, :], grid.zC, color=defaultcolors[j], label=w_lbl, linewidth=lw₁, alpha=α₁, linestyle="-")

    if j > 1
        plot(B_τ[i, :], grid.zC, color=defaultcolors[j], label=τ_lbl, linewidth=lw₂, alpha=α₂, linestyle="--")
        plot(B_b[i, :], grid.zC, color=defaultcolors[j], label=b_lbl, linewidth=lw₂, alpha=α₂, linestyle=":")
    end
    
    sca(axs[2])
    plot(Bz_w[i, :], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    if j > 1
        plot(Bz_τ[i, :], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle="--")
        plot(Bz_b[i, :], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle=":")
    end
    
    sca(axs[3])
    plot(U_w[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    if j > 1
        plot(U_τ[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle="--")
        plot(U_b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle=":")
    end

    sca(axs[4])
    plot(V_w[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    if j > 1
        plot(V_τ[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle="--")
        plot(V_b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle=":")
    end
    
    sca(axs[5])
    plot(w²_w[i, :] / τmax, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    if j > 1
        plot(w²_τ[i, :] / τmax, grid.zF, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle="--")
        plot(w²_b[i, :] / τmax, grid.zF, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, linestyle=":")
    end
end

ylim(-40, 0.1)

sca(axs[1])
xlabel(L"\bar b \, (\mathrm{m \, s^{-2}})")
ylabel("\$ z \$ (m)")

xlim(-0.0004, -0.00005)
legend(loc=3, bbox_to_anchor=(0.05, 0.02, 0.5, 0.5), prop=Dict(:size=>12), frameon=true)

sca(axs[2])
xlabel(L"\partial_z \bar b \, (\mathrm{s^{-2}})")

sca(axs[3])
xlabel(L"\overline{u^\mathrm{L}} \, (\mathrm{m \, s^{-1}})")

sca(axs[4])
xlabel(L"\overline{v^\mathrm{L}} \, (\mathrm{m \, s^{-1}})")

sca(axs[5])
xlabel(L"\overline{\left (w^\mathrm{L} \right )^2} / \max(u_\star^2)")
ylabel("\$ z \$ (m)")

removespines(axs[1], "top", "right")

for ax in axs[2:end-1]
    removespines(ax, "left", "top", "right")
    ax.tick_params(left=false, labelleft=false)
end

removespines(axs[end], "left", "top")
axs[end].tick_params(left=false, labelleft=false, right=true, labelright=true)
axs[end].yaxis.set_label_position("right")

axs[2].set_zorder(-1)
axs[3].set_zorder(-1)

positions = [get_position(ax) for ax in axs]

yshift = 0.05
for (i, pos) in enumerate(positions)
    pos[2] += yshift
    axs[i].set_position(pos)
end

savefig("growth_comparison.png", dpi=480)
