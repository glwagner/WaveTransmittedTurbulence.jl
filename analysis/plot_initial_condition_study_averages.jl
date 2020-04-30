using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2, Printf

ticker = pyimport("matplotlib.ticker")
ScalarFormatter = ticker.ScalarFormatter

fs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

minute = 60.0
hour = 60minute

function plot_initial_condition_study_profiles(wave_amplitude, suffix)
    control_name = "initial_condition_study_resting_0.0x_" * suffix
    resting_name = "initial_condition_study_resting_" * wave_amplitude * "_" * suffix
    excited_name = "initial_condition_study_excited_" * wave_amplitude * "_" * suffix

    control_path = joinpath(@__DIR__, "..", "data", control_name, control_name * "_averages.jld2")
    resting_path = joinpath(@__DIR__, "..", "data", resting_name, resting_name * "_averages.jld2")
    excited_path = joinpath(@__DIR__, "..", "data", excited_name, excited_name * "_averages.jld2")

    grid = get_grid(resting_path)

    control_averages = collect_horizontal_averages(control_path)
    resting_averages = collect_horizontal_averages(resting_path)
    excited_averages = collect_horizontal_averages(excited_path)

    #####
    ##### Plot
    #####

    close("all")
    fig, axs = subplots(ncols=4, figsize=(10, 4), sharey=true)

    α₁ = 0.6
    α₂ = 1.0

    lw₁ = 2.5
    lw₂ = 1.5

    depth = 48

    k = searchsortedfirst(grid.zF, -depth)

    @show length(control_averages.t)
    @show length(resting_averages.t)
    @show length(excited_averages.t)

    iters = [27, 210]

    c_labels, r_labels, e_labels = Tuple([] for i = 1:3)

    f = get_parameter(resting_path, "coriolis", "f")

    for (j, i) in enumerate(iters)
        @show control_averages.t[i] / hour control_averages.t[i] * 1e-4 / 2π

        quarter_time_label = "\$ t = \\frac{1}{4} \\times 2\\pi/f\$"
        full_time_label = "\$ t = 2 \\times 2\\pi/f\$"
        late_time_label = @sprintf("\$ t = %.2f \\times 2\\pi/f\$", f * control_averages.t[i] / 2π)

        if isapprox(f * control_averages.t[i] / 2π, 0.25, atol=1e-2)
            push!(c_labels, "No waves, " * quarter_time_label)
            push!(r_labels, "Resting, " * quarter_time_label)
            push!(e_labels, "Excited, " * quarter_time_label)
        elseif isapprox(f * control_averages.t[i] / 2π, 2, atol=1e-2)
            push!(c_labels, "No waves, " * full_time_label)
            push!(r_labels, "Resting, " * full_time_label)
            push!(e_labels, "Excited, " * full_time_label)
        else
            push!(c_labels, "No waves, " * late_time_label)
            push!(r_labels, "Resting, " * late_time_label)
            push!(e_labels, "Excited, " * late_time_label)
        end
    end

    for (j, i) in enumerate(iters)

        r_label = r_labels[j]
        e_label = e_labels[j]

        sca(axs[1])
        plot(control_averages.b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, ":")  
        plot(resting_averages.b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, "-")
        plot(excited_averages.b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, "--") 

        sca(axs[2])
        plot(control_averages.bz[i, 2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₂, alpha=α₂, ":")  
        plot(resting_averages.bz[i, 2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, "-")
        plot(excited_averages.bz[i, 2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₂, alpha=α₂, "--") 

        sca(axs[3])
        plot(control_averages.S[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, ":")  
        plot(resting_averages.S[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, "-")
        plot(excited_averages.S[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₂, "--") 

        sca(axs[4])
        plot(control_averages.W²[i, :], grid.zF, color=defaultcolors[j], label=c_labels[j], linewidth=lw₂, alpha=α₁, ":")
        plot(resting_averages.W²[i, :], grid.zF, color=defaultcolors[j], label=r_label,     linewidth=lw₁, alpha=α₁, "-")
        plot(excited_averages.W²[i, :], grid.zF, color=defaultcolors[j], label=e_label,     linewidth=lw₂, alpha=α₂, "--") 

    end

    sca(axs[1])
    xlabel(L"\langle b \rangle \, \, (\mathrm{m \, s^{-2}})", labelpad=10.0)
    ylabel("\$ z \$ (m)", labelpad=8.0)

    sca(axs[2])
    xlabel(L"\partial_z \langle b \rangle \, \, (\mathrm{s^{-2}})", labelpad=10.0)

    sca(axs[3])
    xlabel(L"\sqrt{ \left \langle u^\mathrm{L} \right \rangle^2 + \left \langle v^\mathrm{L} \right \rangle^2} \, \, (\mathrm{m \, s^{-1}})")

    sca(axs[4])
    xlabel(L"\left \langle \left ( w^\mathrm{L} \right )^2 \right \rangle \, \, (\mathrm{m^2 \, s^{-2}})")
    legend(loc="lower left", prop=Dict("size" => 12), bbox_to_anchor=(0.27, 0., 1.0, 0.5), frameon=false)

    removespines(axs[1], "top", "right")
    removespines(axs[end], "top", "left")

    for j = 2:4
        removespines(axs[j], "top", "right", "left")
        axs[j].tick_params(left=false, labelleft=false)
    end

    for (ax, lbl) in zip(axs, ("(\\textit{a})", "(\\textit{b})", "(\\textit{c})", "(\\textit{d})"))
        sca(ax)
        text(0.05, 1.025, lbl, transform=ax.transAxes, ha="left", va="bottom")
    end

    positions = [get_position(ax) for ax in axs]

    xwiggle = 0.01
    positions[1][1] -= 2xwiggle
    positions[2][1] -= xwiggle
    positions[3][1] += xwiggle
    positions[4][1] += 2xwiggle

    xshift = 0.015
    yshift = 0.05
    for (i, pos) in enumerate(positions)
        pos[1] -= xshift
        pos[2] += yshift
        axs[i].set_position(pos)
    end

    sca(axs[1])
    xlim(control_averages.b[1, k], control_averages.b[1, end])
    ylim(-depth, 0.1)

    xlim(-7e-4, 0.0)
    axs[1].xaxis.set_ticks([-5e-4, 0.0])

    for ax in axs
        stretch_y!(ax, -0.05)
        stretch_x!(ax, -0.02)
        shift_left!(ax, 0.015)
        shift_up!(ax, 0.02)
    end

    return axs
end
