using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

fs = 16
lblfs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

function name_to_path(name)
    directory = joinpath(@__DIR__, "..", "data", name)
    return joinpath(directory, name * "_averages.jld2"), directory
end

function ∂z(c, grid) 
    cz = zeros(grid.Nz+1)
    cz[2:end-1] = (c[2:end] .- c[1:end-1]) ./ grid.Δz
    return cz
end

function plot_fluxes(suffix, save_points)
     waves_name = "growing_waves_$suffix"
    stress_name = "surface_stress_no_waves_$suffix"
      both_name = "surface_stress_with_waves_$suffix"

     waves_path,  waves_directory = name_to_path(waves_name)
    stress_path, stress_directory = name_to_path(stress_name)
      both_path,   both_directory = name_to_path(both_name)

    grid = get_grid(waves_path)

    wave_amplitude, wavenumber, growth_time_scale = 
        get_multiple_parameters(waves_path, "surface_waves", "wave_amplitude", "wavenumber", "growth_time_scale")

    f = get_parameter(waves_path, "coriolis", "f")

    Uˢ = uˢ(wave_amplitude, wavenumber)
     τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
    u★ = sqrt(τ)

     waves_averages = collect_horizontal_averages(waves_path)
    stress_averages = collect_horizontal_averages(stress_path)
      both_averages = collect_horizontal_averages(both_path)

    for i = 1:length(waves_averages.t)
        @printf "% 5d " i
    end

    @printf "\n"

    for t in waves_averages.t
        @printf "%5.3f " t * f / 2π
    end

    @printf "\n"
    
    #####
    ##### Final profiles
    #####

    close("all")
    fig, axs = subplots(ncols=3, sharey=true, figsize=(10, 5))

    for (j, i) in enumerate(save_points)

         waves_label = @sprintf("waves, \$ t = %.2f f / 2\\pi \$", waves_averages.t[i] * f / 2π)
        stress_label = @sprintf("stress, no waves, \$ t = %.2f f / 2\\pi \$", waves_averages.t[i] * f / 2π)
          both_label = @sprintf("stress with waves, \$ t = %.2f f / 2\\pi \$", waves_averages.t[i] * f / 2π)

        common_kwargs = (color=defaultcolors[j],)

         waves_kwargs = merge((linestyle="-", ), common_kwargs)
        stress_kwargs = merge((linestyle="--",), common_kwargs)
          both_kwargs = merge((linestyle=":", ), common_kwargs)

        @show waves_averages.t[i] * f / 2π

        sca(axs[1])
        plot(∂z( waves_averages.b[:, i], grid), grid.zF,  waves_kwargs..., label=waves_label)
        plot(∂z(stress_averages.b[:, i], grid), grid.zF, stress_kwargs..., label=stress_label)
        plot(∂z(  both_averages.b[:, i], grid), grid.zF,   both_kwargs..., label=both_label)
             
        sca(axs[2])
        plot( waves_averages.wu[:, i] .+  waves_averages.τ₁₃[:, i], grid.zF, waves_kwargs...)
        plot(stress_averages.wu[:, i] .+ stress_averages.τ₁₃[:, i], grid.zF, stress_kwargs...)
        plot(  both_averages.wu[:, i] .+   both_averages.τ₁₃[:, i], grid.zF, both_kwargs...)

        sca(axs[3])
        plot( waves_averages.wv[:, i] .+  waves_averages.τ₂₃[:, i], grid.zF, linestyle="-",  common_kwargs...)
        plot(stress_averages.wv[:, i] .+ stress_averages.τ₂₃[:, i], grid.zF, linestyle="--", common_kwargs...)
        plot(  both_averages.wv[:, i] .+   both_averages.τ₂₃[:, i], grid.zF, linestyle=":",  common_kwargs...)

    end

    return fig, axs
end

#suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"

fig, axs = plot_fluxes(suffix, (27,))# 110, 210,))

sca(axs[1])
legend(loc="lower right")
