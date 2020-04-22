using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Utils

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 16
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

function path_from_name(name)
    directory = joinpath(@__DIR__, "..", "data", name)
         path = joinpath(directory, name * "_averages.jld2")
    return path
end

function plot_hovmoller!(axs, name; velocity_limit=nothing, energy_limit=nothing, stress_limit=nothing,
                                    velocity_saturation=1, energy_saturation=1, stress_saturation=1,
                                    velocity_multiplier=1, energy_multiplier=1, stress_multiplier=1)

    path = path_from_name(name)

    averages = collect_horizontal_averages(path)
    grid = get_grid(path)
    iters = get_iters(path)

    speed = @. sqrt(averages.U^2 + averages.V^2)
    stress = @. averages.wu + averages.τ₁₃
    production = - averages.wu .* averages.Uz - averages.wv .* averages.Vz

    TC, ZC = meshgrid(averages.t / hour, grid.zC)
    TF, ZF = meshgrid(averages.t / hour, grid.zF)

    sca(axs[1])

    velocity_limit, velocity_levels = divergent_limit_levels(averages.U, limit=velocity_limit, 
                                                             saturate=velocity_saturation)

    velocity_contours = contourf(TC, ZC, velocity_multiplier .* averages.U, 
                                 cmap="RdBu", levels=velocity_levels, vmin=-velocity_limit, vmax=velocity_limit)

    sca(axs[2])
    energy_limit, energy_levels = sequential_limit_levels(averages.E, limit=energy_limit,
                                                          saturate=energy_saturation)

    energy_contours = contourf(TC, ZC, energy_multiplier .* averages.E, 
                               vmin=0.0, vmax=energy_limit, levels=energy_levels, cmap="YlGnBu_r")

    sca(axs[3])
    stress_limit, stress_levels = divergent_limit_levels(stress, limit=stress_limit,
                                                         saturate=stress_saturation)

    stress_contours = contourf(TF, ZF, stress_multiplier .* stress, 
                               vmin=-stress_limit, vmax=stress_limit, levels=stress_levels, cmap="RdBu_r")

    return velocity_limit, energy_limit, stress_limit, velocity_contours, energy_contours, stress_contours
end

close("all")
fig, axs = subplots(nrows=3, ncols=2, figsize=(24, 12), sharey=true, sharex=true)

suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"

velocity_limit, energy_limit, stress_limit, velocity_contours, energy_contours, stress_contours =
    plot_hovmoller!(axs[:, 2], "surface_stress_with_waves_" * suffix, velocity_saturation=0.6)

truncate = 10^2
dtick = round(Int, velocity_limit / 3 * truncate) / truncate
velocity_ticks = -3 * dtick:dtick:dtick * 3

truncate = 10^7
dtick = round(Int, stress_limit / 7 * truncate) / truncate
energy_ticks = 0.0 : energy_dtick : energy_limit

truncate = 10^7
dtick = round(Int, stress_limit / 3 * truncate) / truncate
stress_ticks = -3 * dtick:dtick:dtick * 3

velocity_cb = colorbar(velocity_contours, ax=axs[1, 2], shrink=0.8, aspect=10, pad=0.05, 
                       ticks=velocity_ticks)

energy_cb = colorbar(energy_contours, ax=axs[2, 2], shrink=0.8, aspect=10, pad=0.05, format="%.2e")
                     #ticks=energy_ticks, 

stress_cb = colorbar(stress_contours, ax=axs[3, 2], shrink=0.8, aspect=10, pad=0.05, 
                     ticks=stress_ticks, format="%.2e")

velocity_limit, energy_limit, stress_limit, _, energy_contours, stress_contours =
    plot_hovmoller!(axs[:, 1], "growing_waves_" * suffix, velocity_limit=velocity_limit)

#truncate = 10^2
#dtick = round(Int, velocity_limit / 3 * truncate) / truncate
#velocity_ticks = -3 * dtick:dtick:dtick * 3

truncate = 10^6
dtick = round(Int, energy_limit / 7 * truncate) / truncate
energy_ticks = 0.0 : energy_dtick : energy_limit

truncate = 10^9
dtick = round(Int, stress_limit / 3 * truncate) / truncate
stress_ticks = -3 * dtick:dtick:dtick * 3

velocity_cb = colorbar(velocity_contours, ax=[axs[1, 1]], shrink=0.8, aspect=10, ticks=velocity_ticks,
                       pad=0.05, location="left")

energy_cb = colorbar(energy_contours, ax=[axs[2, 1]], shrink=0.8, aspect=10,# ticks=energy_ticks,
                     pad=0.05, location="left", format="%.2e")

stress_cb = colorbar(stress_contours, ax=[axs[3, 1]], shrink=0.8, aspect=10, ticks=stress_ticks,
                     pad=0.05, location="left", format="%.2e")

ylim(-48, 0)
xlim(0, 24)

for j = 1:2
    bottom_left_text!(axs[1, j], "\$ U \$",                                             color="k", fontsize=18)
    bottom_left_text!(axs[2, j], "\$ E \$",                                             color="w", fontsize=18)
    bottom_left_text!(axs[3, j], "\$ - \\langle w' u' \\rangle + \\mathcal{T}_{xz} \$", color="k", fontsize=18)
end

sca(axs[1, 1])
xlabel("time (hours)"; labelpad=12.0)
axs[1, 1].xaxis.set_label_position("top")

sca(axs[1, 2])
xlabel("time (hours)"; labelpad=12.0)
axs[1, 2].xaxis.set_label_position("top")

sca(axs[3, 1])
xlabel("time (hours)")

sca(axs[3, 2])
xlabel("time (hours)")

axs[1, 1].tick_params(left=false, labelleft=false, right=true, labelright=true, top=true, labeltop=true)
axs[1, 2].tick_params(left=true,  labelleft=true,  top=true, labeltop=true)
axs[2, 1].tick_params(left=false, labelleft=false, right=true, labelright=true)
axs[3, 1].tick_params(left=false, labelleft=false, right=true, labelright=true)

axs[2, 2].tick_params(left=true,  labelleft=true)
axs[3, 2].tick_params(left=true,  labelleft=true)

for i = 1:3
    sca(axs[i, 1])
    axs[i, 1].yaxis.set_label_position("right")
    ylabel(" \$ z \$ (meters)", labelpad=28.0)
end

savefig(joinpath(@__DIR__, "..", "figures", "compare_hovmoller.png"), dpi=480)
