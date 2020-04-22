using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Utils

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 12
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

name = "surface_stress_with_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
#name = "growing_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"

directory = joinpath(@__DIR__, "..", "data", name)
     path = joinpath(directory, name * "_averages.jld2")

file = jldopen(path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

averages = collect_horizontal_averages(path)
grid = get_grid(path)
iters = get_iters(path)

speed = @. sqrt(averages.U^2 + averages.V^2)
stress = @. averages.wu + averages.τ₁₃
production = - averages.wu .* averages.Uz - averages.wv .* averages.Vz

TC, ZC = meshgrid(averages.t / hour, grid.zC)
TF, ZF = meshgrid(averages.t / hour, grid.zF)

function sequential_limit_levels(field, saturate=0.8; nlevs=15)
    max_field = maximum(field)
    limit = saturate * max_field
    levels = vcat(range(0.0, stop=limit, length=nlevs), [max_field])
    return limit, levels
end

function divergent_limit_levels(field, saturate=0.6; nlevs=14)
    max_field = maximum(abs, field)
    limit = saturate * max_field
    levels = vcat([-max_field], range(-limit, stop=limit, length=nlevs), [max_field])
    return limit, levels
end

function bottom_left_text!(ax, txt; kwargs...)
    sca(ax)
    text(0.01, 0.02, txt; transform=ax.transAxes, horizontalalignment="left", verticalalignment="bottom", kwargs...)
    return nothing
end

close("all")
fig, axs = subplots(nrows=4, ncols=2, figsize=(24, 12), sharey=true, sharex=true)

sca(axs[1, 1])
energy_limit, energy_levels = sequential_limit_levels(averages.E)
contourf(TC, ZC, averages.E, vmin=0.0, vmax=energy_limit, levels=energy_levels, cmap="YlGnBu_r")

axs[1, 1].tick_params(labeltop=true, top=true)
axs[1, 2].tick_params(labeltop=true, top=true)

sca(axs[2, 1])
production_limit, production_levels = divergent_limit_levels(production)
contourf(TF, ZF, production, vmin=-production_limit, vmax=production_limit, levels=production_levels, cmap="RdBu_r")

sca(axs[3, 1])
stress_limit, stress_levels = divergent_limit_levels(stress)
contourf(TF, ZF, stress, vmin=-stress_limit, vmax=stress_limit, levels=stress_levels, cmap="RdBu_r")

sca(axs[4, 1])
variance_limit, variance_levels = sequential_limit_levels(averages.W²)
contourf(TF, ZF, averages.W², cmap="YlGnBu_r", levels=variance_levels, vmin=0.0, vmax=variance_limit)

#####
##### Second column
#####

sca(axs[1, 2])
mean_energy_limit, mean_energy_levels = sequential_limit_levels(speed.^2)
contourf(TC, ZC, speed.^2, vmin=0.0, vmax=mean_energy_limit, levels=mean_energy_levels, cmap="YlGnBu_r")

sca(axs[2, 2])
velocity_limit, velocity_levels = divergent_limit_levels(averages.U)
contourf(TC, ZC, averages.U, cmap="RdBu", levels=velocity_levels, vmin=-velocity_limit, vmax=velocity_limit)

sca(axs[3, 2])
shear_limit, shear_levels = divergent_limit_levels(averages.Uz)
contourf(TF, ZF, averages.Uz, cmap="RdBu", levels=shear_levels, vmin=-shear_limit, vmax=shear_limit)

sca(axs[4, 2])
triple_limit, triple_levels = divergent_limit_levels(averages.W³)
contourf(TF, ZF, averages.W³, cmap="RdBu_r", levels=triple_levels, vmin=-triple_limit, vmax=triple_limit)

ylim(-48, 0)
xlim(0, 24)

bottom_left_text!(axs[1, 1], "\$ E \$",                                             color="w", fontsize=18)
bottom_left_text!(axs[2, 1], "\$ - \\langle w' u'_i \\rangle \\partial_z U_i \$",   color="k", fontsize=18)
bottom_left_text!(axs[3, 1], "\$ - \\langle w' u' \\rangle + \\mathcal{T}_{xz} \$", color="k", fontsize=18)
bottom_left_text!(axs[4, 1], "\$ \\langle w'^2 \\rangle \$",                        color="w", fontsize=18)

bottom_left_text!(axs[1, 2], "\$ \\frac{1}{2} \\left ( U^2 + V^2 \\right ) \$", color="w", fontsize=18)
bottom_left_text!(axs[2, 2], "\$ U \$",                                         color="k", fontsize=18)
bottom_left_text!(axs[3, 2], "\$ \\partial_z U \$",                             color="k", fontsize=18)
bottom_left_text!(axs[4, 2], "\$ \\langle w'^3 \\rangle \$",                    color="k", fontsize=18)

sca(axs[1, 1])
xlabel("time (hours)"; labelpad=12.0)
axs[1, 1].xaxis.set_label_position("top")

sca(axs[1, 2])
xlabel("time (hours)"; labelpad=12.0)
axs[1, 2].xaxis.set_label_position("top")

sca(axs[4, 1])
xlabel("time (hours)")

sca(axs[4, 2])
xlabel("time (hours)")

for i = 1:4
    sca(axs[i, 1])
    ylabel(" \$ z \$ (meters)")
end
