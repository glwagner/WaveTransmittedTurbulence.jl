using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

     run_name = "growing_wave_forced_Qb1.0e-09_a1.5_k6.3e-02_T4.0_Nh256_Nz256"
run_directory = joinpath(@__DIR__, "..", "data", run_name)
     run_path = joinpath(run_directory, run_name * "_fields_part3.jld2")

file = jldopen(run_path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
      wave_number = file["surface_waves/wave_number"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

Uˢ = uˢ(wave_amplitude, wave_number)
 τ = wave_amplitude^2 * √(g_Earth * wave_number) / (2 * growth_time_scale)
u★ = sqrt(τ)
grid = get_grid(run_path)

iters = get_iters(run_path)

t, U_w, V_w, B_w, Bz_w, w²_w, E_w = extract_averages_timeseries(run_directory)

#####
##### Final profiles
#####

fig, axs = subplots(ncols=5, sharey=true, figsize=(10, 4))

lw₁ = 2
lw₂ = 1.5
α₁ = 0.6
α₂ = 1.0

i = 1
j = 1
w_lbl = "Initial_condition"

sca(axs[1])
plot(B_w[i][2:end-1], grid.zC, color=defaultcolors[j], label=w_lbl, linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[2])
plot(Bz_w[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[3])
plot(U_w[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[4])
plot(V_w[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[5])
plot(w²_w[i][2:end-1] / τ, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

i = 5 #length(B_w)
j = 2

sca(axs[1])
plot(B_w[i][2:end-1], grid.zC, color=defaultcolors[j], label=w_lbl, linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[2])
plot(Bz_w[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[3])
plot(U_w[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[4])
plot(V_w[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[5])
plot(w²_w[i][2:end-1] / τ, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")



#=
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
=#
