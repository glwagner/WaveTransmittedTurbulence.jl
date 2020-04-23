using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Utils

using Oceananigans.Buoyancy: g_Earth

using Statistics

fs = 16
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

function path_from_name(name)
    directory = joinpath(@__DIR__, "..", "data", name)
         path = joinpath(directory, name * "_averages.jld2")
    return path
end

integrate_profile_timeseries(profile, grid) = dropdims(mean(profile, dims=2), dims=2)

function plot_volume_budget!(axs, name; kwargs...)

    path = path_from_name(name)

    averages = collect_horizontal_averages(path)
    grid = get_grid(path)
    iters = get_iters(path)

    mean_energy_timeseries = @. (averages.U^2 + averages.V^2) / 2

    production_faces = - averages.wu .* averages.Uz - averages.wv .* averages.Vz
    production = @. (production_faces[:, 2:end] + production_faces[:, 1:end-1]) / 2

    ℰ = integrate_profile_timeseries(mean_energy_timeseries, grid)
    E = integrate_profile_timeseries(averages.E, grid)
    υ = integrate_profile_timeseries(averages.W², grid) ./ 2
    P = integrate_profile_timeseries(production, grid)

    sca(axs[2])
    U = integrate_profile_timeseries(averages.U, grid)
    V = integrate_profile_timeseries(averages.V, grid)
    plot(averages.t / hour, U, "k--", alpha=0.6, linewidth=1.5)
    plot(averages.t / hour, V, "k:", alpha=0.6, linewidth=1.5)
    xticks(0:4:32)

    sca(axs[3])
    plot(averages.t / hour, ℰ .+ E,  "-"; kwargs...)
    xticks(0:4:32)

    sca(axs[4])
    plot(averages.t / hour, E, "--"; kwargs...)
    xticks(0:4:32)

    sca(axs[5])
    plot(averages.t / hour, P, "-"; kwargs...)
    xticks(0:4:32)

    return nothing
end

suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
name = "growing_waves_" * suffix
path = path_from_name(name)

file = jldopen(path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]

close(file)

averages = collect_horizontal_averages(path)
grid = get_grid(path)
U = integrate_profile_timeseries(averages.U, grid)
V = integrate_profile_timeseries(averages.V, grid)

# Effective stress
τʷ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
τ(t) = τʷ * t / growth_time_scale * exp(-t^2 / (2*growth_time_scale^2))

close("all")
fig, axs = subplots(nrows=5, sharex=true, figsize=(8, 10))

suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"

sca(axs[1])
plot(averages.t / hour, τ.(averages.t), "-")
xticks(0:4:32)

plot_volume_budget!(axs, "growing_waves_" * suffix,             label="Growing waves")
plot_volume_budget!(axs, "surface_stress_with_waves_" * suffix, label="Surface stress, steady waves")
plot_volume_budget!(axs, "surface_stress_no_waves_" * suffix,   label="Surface stress, no waves")

removespines(axs[1], "right")
removespines(axs[2], "top", "right")
removespines(axs[3], "top", "right")
removespines(axs[4], "top", "right")
removespines(axs[5], "top", "right")

axs[1].tick_params(top=true, labeltop=true)

sca(axs[1])
ylabel(L"\int \partial_t u^\mathrm{S} \, \mathrm{d} z \, \, \mathrm{(m^2 \, s^{-2})}")
xlabel("time (hours)", labelpad=12.0)
axs[1].xaxis.set_label_position("top")

sca(axs[2])
ylabel(L"\left \langle u_i \right \rangle \,  \, \mathrm{(m \, s^{-1})}")

sca(axs[3])
ylabel(L"\left \langle \frac{1}{2} u_i^2 \right \rangle \,  \, \mathrm{(m^2 \, s^{-2})}")

sca(axs[4])
ylabel(L"\left \langle \frac{1}{2} u_i'^2 \right \rangle \,  \, \mathrm{(m^2 \, s^{-2})}")
legend(bbox_to_anchor=(0.5, 0.25, 1.0, 1.0), loc=2, frameon=false)

sca(axs[5])
ylabel(L"-H^{-1} \int \overline{w'u_i'} \partial_z U_i \, \mathrm{d} z") # \,  \, \mathrm{(m^2 \, s^{-3})}")
xlabel("time (hours)")

xlim(-1, 32)
#tight_layout()

function shift_right!(ax, shift=0.01)
    pos = get_position(ax)
    pos[1] += shift
    ax.set_position(pos)
    return nothing
end

function shift_up!(ax, shift=0.01)
    pos = get_position(ax)
    pos[2] += shift
    ax.set_position(pos)
    return nothing
end

for (i, ax) in enumerate(axs)
    shift_right!(ax, 0.05)
    shift_up!(ax, 0.05)
end

for (i, ax) in enumerate(axs)
    shift_up!(ax, -0.02 * (i-1))
end

savefig(joinpath(@__DIR__, "..", "figures", "budget_timeseries.png"), dpi=480)
