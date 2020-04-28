using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2, Printf

function plot_initial_condition_study_profiles(wave_amplitude, suffix)
    control_name = "initial_condition_study_resting_0.0x_" * suffix
    resting_name = "initial_condition_study_resting_" * wave_amplitude * "_" * suffix
    resting_name = "initial_condition_study_excited_" * wave_amplitude * "_" * suffix

    control_path = joinpath(@__DIR__, "..", "data", control_name, control_name * "_fields_part4.jld2")
    resting_path = joinpath(@__DIR__, "..", "data", resting_name, resting_name * "_fields_part4.jld2")
    excited_path = joinpath(@__DIR__, "..", "data", excited_name, excited_name * "_fields_part4.jld2")

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

    lw₂ = 2.5
    lw₁ = 1.5

    depth = 36

    k = searchsortedfirst(grid.zF, -depth)

    @show length(control_averages.t)
    @show length(resting_averages.t)
    @show length(excited_averages.t)

    iters = [51, 401]

    c_labels, r_labels, e_labels = Tuple([] for i = 1:3)

    f = get_parameter(resting_path, "coriolis", "f")

    for (j, i) in enumerate(iters)
        @show control.t[i] / hour control.t[i] * 1e-4 / 2π

        if j == 1
            if isapprox(f * control.t[i] / 2π, 0.25, atol=1e-2)
                push!(c_labels, @sprintf("No waves, \$ t = \\frac{1}{4} \\times 2\\pi/f\$"))
            else
                push!(c_labels, @sprintf("No waves, \$ t = %.2f \\times 2\\pi/f\$", f * control.t[i] / 2π))
            end
        else
            if isapprox(f * control.t[i] / 2π, 2.0, atol=1e-2)
                push!(c_labels, @sprintf("\$ t = %d \\times 2\\pi/f\$", f * control.t[i] / 2π))
            else
                push!(c_labels, @sprintf("\$ t = %.2f \\times 2\\pi/f\$", f * control.t[i] / 2π))
            end
        end

        push!(r_labels, "Equilibrated")
        push!(e_labels, "Excited")
    end

    for (j, i) in enumerate(iters)

        r_label = r_labels[j]
        e_label = e_labels[j]

        sca(axs[1])
        plot(control_averages.b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₁, "-")  
        plot(resting_averages.b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, "--")
        plot(excited_averages.b[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, ":") 

        sca(axs[2])
        plot(control_averages.bz[i, 2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₂, alpha=α₁, "-")
        plot(resting_averages.bz[i, 2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₂, "--")
        plot(excited_averages.bz[i, 2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₂, ":")

        sca(axs[3])
        plot(control_averages.U[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₁, label=c_labels[j], "-")
        plot(resting_averages.U[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, label=r_label, "--")
        plot(excited_averages.U[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, label=e_label, ":")

        sca(axs[4])
        plot(control_averages.V[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₁, "-")
        plot(resting_averages.V[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, "--")
        plot(excited_averages.V[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, ":")

    end

    sca(axs[1])
    xlabel(L"\bar b \, (\mathrm{m \, s^{-2}})")
    ylabel("\$ z \$ (m)")

    sca(axs[2])
    xlabel(L"\partial_z \bar b \, (\mathrm{s^{-2}})")

    sca(axs[3])
    xlabel(L"\bar u \, (\mathrm{m \, s^{-1}})")
    legend(loc=3, prop=Dict("size" => 12), bbox_to_anchor=(0.3, 0., 1.0, 0.5))

    sca(axs[4])
    xlabel(L"\bar v \, (\mathrm{m \, s^{-1}})")
    ylabel("\$ z \$ (m)")

    removespines(axs[1], "top", "right")
    removespines(axs[end], "top", "left")

    for j = 2:3
        removespines(axs[j], "top", "right", "left")
        axs[j].tick_params(left=false, labelleft=false)
    end

    axs[end].yaxis.set_label_position("right")
    axs[end].tick_params(left=false, labelleft=false, right=true, labelright=true)
    axs[end].set_zorder(-1)

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

    #axs[1].xaxis.set_ticks([-5e-4, -2e-4]) #set_major_formatter(formatter)

    formatter = ScalarFormatter() #useMathText=true)
    formatter.set_powerlimits((-2, 2))
    axs[1].xaxis.set_major_formatter(formatter)

    formatter = ScalarFormatter() #useMathText=true)
    formatter.set_powerlimits((-2, 2))
    axs[2].xaxis.set_major_formatter(formatter)

    return axs
end

axs = plot_initial_condition_study_profiles("1.0x", "Nh256_Nz256")
savefig(joinpath(@__DIR__, "..", "figures", "figure_7.png"), dpi=480)

axs = plot_initial_condition_study_profiles("4.0x", "Nh256_Nz256")
savefig(joinpath(@__DIR__, "..", "figures", "figure_8.png"), dpi=480)
