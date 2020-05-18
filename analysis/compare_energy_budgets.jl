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

function plot_volume_budget!(axs, name; ∂t_uˢ=nothing, τ=nothing, kwargs...)

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

    U = integrate_profile_timeseries(averages.U, grid)
    V = integrate_profile_timeseries(averages.V, grid)

    power = nothing

    if τ != nothing
        power = @. τ(averages.t) * averages.U[:, end]
    elseif ∂t_uˢ != nothing
        t = reshape(averages.t, length(averages.t), 1)
        zC = reshape(grid.zC, 1, grid.Nz)
        depth_distributed_power = @. ∂t_uˢ(zC, t) * averages.U
        power = grid.Lz .* integrate_profile_timeseries(depth_distributed_power, grid)
    else
        power = zeros(length(averages.t))
    end

    sca(axs[2])
    plot(averages.t / hour, U, "k--", alpha=0.6, linewidth=1.5)
    plot(averages.t / hour, V, "k:", alpha=0.6, linewidth=1.5)
    xticks(0:4:32)

    sca(axs[3])
    plot(averages.t / hour, power, "-"; kwargs...)
    xticks(0:4:32)

    sca(axs[4])
    plot(averages.t / hour, ℰ,  "-"; kwargs...)
    xticks(0:4:32)

    sca(axs[5])
    plot(averages.t / hour, E, "--"; kwargs...)
    xticks(0:4:32)

    sca(axs[6])
    plot(averages.t / hour, P, "-"; kwargs...)
    xticks(0:4:32)

    return nothing
end

suffix = "Qb5.0e-10_Nsq1.0e-06_f1.0e-04_dom1.5_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz384"
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
∂t_uˢ(z, t) = exp(2 * wavenumber * z) * wave_amplitude^2 * wavenumber * √(g_Earth * wavenumber) * t / growth_time_scale^2 * exp(-t^2 / (2*growth_time_scale^2))

close("all")
fig, axs = subplots(nrows=6, sharex=true, figsize=(8, 10))

sca(axs[1])
plot(averages.t / hour, τ.(averages.t), "-")
xticks(0:4:32)

plot_volume_budget!(axs, "growing_waves_" * suffix,                    ∂t_uˢ=∂t_uˢ, color=defaultcolors[1], label="Growing waves")
plot_volume_budget!(axs, "surface_stress_with_steady_waves_" * suffix, τ=τ,         color=defaultcolors[2], label="Surface stress, steady waves")
plot_volume_budget!(axs, "surface_stress_no_waves_" * suffix,          τ=τ,         color=defaultcolors[3], label="Surface stress, no waves")

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
ylabel("power")

sca(axs[4])
ylabel(L"\frac{1}{2} \left \langle u_i \right \rangle^2 \,  \, \mathrm{(m^2 \, s^{-2})}")

sca(axs[5])
ylabel(L"\left \langle \frac{1}{2} u_i'^2 \right \rangle \,  \, \mathrm{(m^2 \, s^{-2})}")
legend(bbox_to_anchor=(0.5, 0.25, 1.0, 1.0), loc=2, frameon=false)

sca(axs[6])
ylabel(L"-\left \langle \overline{w'u_i'} \partial_z U_i \right \rangle") # \,  \, \mathrm{(m^2 \, s^{-3})}")
xlabel("time (hours)")

xlim(-1, 32)

for (i, ax) in enumerate(axs)
    shift_right!(ax, 0.05)
    shift_up!(ax, 0.05)
end

for (i, ax) in enumerate(axs)
    shift_up!(ax, -0.02 * (i-1))
end

savefig(joinpath(@__DIR__, "..", "figures", "budget_timeseries.png"), dpi=480)
