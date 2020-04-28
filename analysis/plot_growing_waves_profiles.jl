using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

fs = 13
lblfs = 13
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

function name_to_path(name)
    directory = joinpath(@__DIR__, "..", "data", name)
    return joinpath(directory, name * "_averages.jld2"), directory
end

function plot_growing_wave_profiles(suffix; i=4200, ylims=nothing) #(-64, 0.1))
    no_stress_growing_waves_name = "growing_waves_$suffix"
    stress_no_waves_name = "surface_stress_no_waves_$suffix"
    #stress_steady_waves_name = "surface_stress_with_steady_waves_$suffix"
    stress_steady_waves_name = "surface_stress_with_waves_$suffix"
    stress_growing_waves_name = "surface_stress_with_growing_waves_$suffix"

    no_stress_growing_waves_path,  no_stress_growing_waves_directory = 
        name_to_path(no_stress_growing_waves_name)

    stress_no_waves_path,  stress_no_waves_directory = 
        name_to_path(stress_no_waves_name)

    stress_steady_waves_path,  stress_steady_waves_directory = 
        name_to_path(stress_steady_waves_name)

    stress_growing_waves_path,  stress_growing_waves_directory = 
        name_to_path(stress_growing_waves_name)

    grid = get_grid(no_stress_growing_waves_path)

    wave_amplitude, wavenumber, growth_time_scale = 
        get_multiple_parameters(no_stress_growing_waves_path, 
                                "surface_waves", "wave_amplitude", "wavenumber", "growth_time_scale")

    f = get_parameter(no_stress_growing_waves_path, "coriolis", "f")

    Uˢ = uˢ(wave_amplitude, wavenumber)
     τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
    u★ = sqrt(τ)

    no_stress_growing_waves = collect_horizontal_averages(no_stress_growing_waves_path)
    stress_no_waves = collect_horizontal_averages(stress_no_waves_path)
    stress_steady_waves = collect_horizontal_averages(stress_steady_waves_path)
    stress_growing_waves = collect_horizontal_averages(stress_growing_waves_path)

    @show dt = (no_stress_growing_waves.t[2]-no_stress_growing_waves.t[1]) * f/2π
    @show 1/dt

    @show length(no_stress_growing_waves.t)
    @show no_stress_growing_waves.t[i] * f/2π 

    #####
    ##### Final profiles
    #####

    close("all")
    fig, axs = subplots(ncols=4, sharey=true, figsize=(10.2, 4.2))

    plot_profiles!(axs, grid, 
                   no_stress_growing_waves.b[1, :], 
                   no_stress_growing_waves.bz[1, :], 
                   nothing, 
                   no_stress_growing_waves.W²[1, :] / τ, 
                   label="Initial condition", color="k", linewidth=1.5, alpha=0.4, linestyle="--")

    plot_profiles!(axs, grid, 
                   no_stress_growing_waves.b[i, :],
                   no_stress_growing_waves.bz[i, :], 
                   no_stress_growing_waves.S[i, :], 
                   no_stress_growing_waves.W²[i, :] / τ, 
                   label="Growing swell, no surface stress, \$ t = 2 \\pi / f \$", 
                   linestyle="-", color=defaultcolors[1], alpha=1.0, linewidth=2)

    plot_profiles!(axs, grid, 
                   stress_no_waves.b[i, :],
                   stress_no_waves.bz[i, :], 
                   stress_no_waves.S[i, :], 
                   stress_no_waves.W²[i, :] / τ, 
                   label="Surface stress, no swell, \$ t = 2 \\pi / f \$", 
                   linestyle="-", color=defaultcolors[2], alpha=0.8, linewidth=1.5)

    plot_profiles!(axs, grid, 
                   stress_steady_waves.b[i, :],
                   stress_steady_waves.bz[i, :], 
                   stress_steady_waves.S[i, :], 
                   stress_steady_waves.W²[i, :] / τ, 
                   label="Surface stress, steady swell, \$ t = 2 \\pi / f \$", 
                   linestyle="-", color="xkcd:red", alpha=0.8, linewidth=1.5)

    plot_profiles!(axs, grid, 
                   stress_growing_waves.b[i, :],
                   stress_growing_waves.bz[i, :], 
                   stress_growing_waves.S[i, :], 
                   stress_growing_waves.W²[i, :] / τ, 
                   label="Surface stress, growing swell, \$ t = 2 \\pi / f \$", 
                   linestyle="-", color="xkcd:aqua green", alpha=0.8, linewidth=1.5)

    if ylims != nothing
        ylim(ylims...)
    end

    sca(axs[1])
    xlabel(L"\langle b \rangle \, \, \, (\mathrm{m \, s^{-2}})", fontsize=lblfs, labelpad=10.0)
    ylabel("\$ z \$ (m)", fontsize=lblfs, labelpad=8.0)

    xlim(-6e-5, 0.0)
    xticks([-5e-5, -1e-5], [L"-5 \times 10^{-5}", L"-10^{-5}"])
    legend(loc="lower left", bbox_to_anchor=(0.01, 0.01, 0.99, 0.99), prop=Dict(:size=>lblfs), ncol=2, frameon=true, framealpha=0.9)

    sca(axs[2])

    xlim(-1e-6, 4e-6)

    xticks([0, 2e-6], ["0", L"2 \times 10^{-6}"])
    xlabel(L"\partial_z \langle b \rangle \, \, \, (\mathrm{s^{-2}})", fontsize=lblfs, labelpad=10.0)

    sca(axs[3])
    xlabel(L"\sqrt{\left \langle u^\mathrm{L} \right \rangle^2 + \left \langle v^\mathrm{L} \right \rangle^2} \, \, \, (\mathrm{m \, s^{-1}})", 
           fontsize=lblfs, labelpad=4.0)

    sca(axs[4])
    xlabel(L"\left \langle w^\mathrm{L} \right \rangle^2 \, / \, \max(u_\star^2)", fontsize=lblfs, labelpad=8.0)
    ylabel("\$ z \$ (m)", fontsize=lblfs, labelpad=10.0)

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
    axs[4].set_zorder(-1)

    for (ax, lbl) in zip(axs, ("(\\textit{a})", "(\\textit{b})", "(\\textit{c})", "(\\textit{d})"))
        shift_up!(ax, 0.05)
        shift_left!(ax, 0.02)
        sca(ax)
        text(0.05, 1.02, lbl, transform=ax.transAxes, fontsize=lblfs, ha="left", va="bottom")
    end

    return fig, axs
end
