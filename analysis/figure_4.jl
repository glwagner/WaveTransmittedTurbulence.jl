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

function plot_depth_time!(axs, name; velocity_limit=nothing, energy_limit=nothing, stress_limit=nothing,
                                    velocity_saturation=1, energy_saturation=1,
                                    velocity_multiplier=1, energy_multiplier=1)

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
                                 cmap="RdBu_r", levels=velocity_levels, vmin=-velocity_limit, vmax=velocity_limit)

    sca(axs[2])
    energy_limit, energy_levels = sequential_limit_levels(averages.E, limit=energy_limit,
                                                          saturate=energy_saturation)

    energy_contours = contourf(TC, ZC, energy_multiplier .* averages.E, 
                               vmin=0.0, vmax=energy_limit, levels=energy_levels, cmap="YlGnBu_r")

    return velocity_limit, energy_limit, velocity_contours, energy_contours
end

close("all")
fig, axs = subplots(nrows=3, ncols=2, figsize=(15, 9), sharex=true)

suffix = "Qb5.0e-10_Nsq1.0e-06_f1.0e-04_dom1.5_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz384"

#####
##### Right side: surface stress with steady waves
#####

velocity_limit, energy_limit, velocity_contours, energy_contours =
    plot_depth_time!(axs[2:3, 2], "surface_stress_with_steady_waves_" * suffix, velocity_saturation=0.6, energy_saturation=0.5)

dtick = velocity_limit / 2
velocity_ticks = -2 * dtick : dtick : dtick * 2

dtick = energy_limit / 4
energy_ticks = 0.0 : dtick : energy_limit

shrink, aspect, pad = 0.8, 8, 0.0
velocity_cb = colorbar(velocity_contours, ax=axs[2, 2], shrink=shrink, aspect=aspect, pad=pad, ticks=velocity_ticks, format="%.2f")
  energy_cb = colorbar(energy_contours,   ax=axs[3, 2], shrink=shrink, aspect=aspect, pad=pad, ticks=energy_ticks, format="%.1e")

velocity_cb.ax.set_title(L"U^\mathrm{L} \, \, \mathrm{(m \, s^{-1})}", fontsize=16, pad=12.0)
energy_cb.ax.set_title(L"E \, \, \mathrm{(m^2 \, s^{-2})}",            fontsize=16, pad=12.0)

#####
##### Left side: growing waves
#####

velocity_limit, energy_limit, _, energy_contours =
    plot_depth_time!(axs[2:3, 1], "growing_waves_" * suffix, velocity_limit=velocity_limit, energy_limit=energy_limit)

bottom_left_text!(axs[2, 1], "Growing waves", color="k", fontsize=16)
bottom_left_text!(axs[2, 2], "Surface stress with steady waves", color="k", fontsize=16)

bottom_left_text!(axs[3, 1], "Growing waves", color="w", fontsize=16)
bottom_left_text!(axs[3, 2], "Surface stress with steady waves", color="w", fontsize=16)

sca(axs[3, 1])
xlabel("time (hours)")

sca(axs[3, 2])
xlabel("time (hours)")

for ax in axs
    ax.tick_params(left=true, right=true, bottom=true, top=true)
end

axs[1, 1].tick_params(top=false, labelbottom=true, right=false)
axs[1, 2].tick_params(top=false, labelbottom=true, left=false, labelleft=false, labelright=true)

for i = (2, 3)
    for j = (1, 2)
        sca(axs[i, j])
        ylabel(" \$ z \$ (meters)", labelpad=12.0)
    end
    axs[i, 2].yaxis.set_label_position("right")
    axs[i, 2].tick_params(labelleft=false, labelright=true)
end


#####
##### Plot effective stress
#####

path = path_from_name("growing_waves_" * suffix)
averages = collect_horizontal_averages(path)
file = jldopen(path)
   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]
close(file)

τw = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
τ(t) = τw * t / growth_time_scale * exp(-t^2 / (2 * growth_time_scale^2))

function plot_effective_stress()
    plot(averages.t / hour, τ.(averages.t))

    plot([1.0, 1.0] * π/f / hour,  [2.15e-5, 4.1e-5], "k", alpha=0.4, linewidth=1, zorder=1)
    plot([1.0, 1.0] * 2π/f / hour, [1e-6,  2e-5], "k", alpha=0.4, linewidth=1, zorder=2)

    text(π/f / hour,  5e-5, L"\pi/f", fontsize=fs-2,    ha="center", va="bottom")
    text(2π/f / hour, 3e-5, L"2 \pi/f", fontsize=fs-2, ha="center", va="bottom")

    xlabel("time (hours)", labelpad=8.0)
    yticks([0, 5e-5], ["0", L"5 \times 10^{-5}"])
    return nothing
end

sca(axs[1, 1])
plot_effective_stress()
removespines("top", "right")
ylabel(L"\partial_t U^\mathrm{S} \, \, (\mathrm{m^2 \, s^{-2}})", labelpad=12.0)

sca(axs[1, 2])
plot_effective_stress()
removespines("top", "left")
axs[1, 2].yaxis.set_label_position("right")
ylabel(L"\tau \, \, (\mathrm{m^2 \, s^{-2}})", labelpad=12.0)

pos2 = get_position(axs[2, 1])
pos1 = get_position(axs[1, 1])
pos1[1] = pos2[1]
pos1[3] = pos2[3]
pos1[2] += 0.07
pos1[4] /= 2
axs[1, 1].set_position(pos1)

pos2 = get_position(axs[2, 2])
pos1 = get_position(axs[1, 2])
pos1[1] = pos2[1]
pos1[3] = pos2[3]
pos1[2] += 0.07
pos1[4] /= 2
axs[1, 2].set_position(pos1)

for ax in axs
    sca(ax)
    xlim(0, 24)
    xticks([0, 4, 8, 12, 16, 20, 24])
end

for i = 2:3
    sca(axs[i, 1])
    ylim(-48, 0)

    sca(axs[i, 2])
    ylim(-48, 0)
end

for i = 1:3
    pos1 = get_position(axs[i, 1])
    pos2 = get_position(axs[i, 2])
    Δx = pos1[3] - pos2[3]
    frac = 1/2
    stretch_x!(axs[i, 2], frac * Δx)
    stretch_x!(axs[i, 1], (frac - 1) * Δx)
end

for i = 1:3
    shift_left!(axs[i, 1], 0.02)
    shift_left!(axs[i, 2], 0.06)
end

for ax in axs
    shift_up!(ax, 0.05)
end

for cb in (velocity_cb, energy_cb)
    for tick in cb.ax.get_yticklabels()
        tick.set_fontsize(14)
    end
    shift_right!(cb.ax, 0.065)
    shift_up!(cb.ax, 0.04)
end

savefig(joinpath(@__DIR__, "..", "figures", "figure_4.png"), dpi=480)
