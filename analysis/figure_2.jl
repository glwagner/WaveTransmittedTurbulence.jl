using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

suffix = "Qb1.0e-09_a2.0_k6.3e-02_T4.0_Nh256_Nz256" 
part = 2

   waves_name = "growing_waves_$suffix"
  stress_name = "surface_stress_no_waves_$suffix"
    both_name = "surface_stress_with_waves_$suffix"

 waves_directory = joinpath(@__DIR__, "..", "data", waves_name)
stress_directory = joinpath(@__DIR__, "..", "data", stress_name)
  both_directory = joinpath(@__DIR__, "..", "data", both_name)

 waves_path = joinpath(waves_directory,   waves_name * "_fields_part$part.jld2")
stress_path = joinpath(stress_directory, stress_name * "_fields_part$part.jld2")
  both_path = joinpath(both_directory,     both_name * "_fields_part$part.jld2")

file = jldopen(waves_path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

Uˢ = uˢ(wave_amplitude, wavenumber)
 τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
u★ = sqrt(τ)
grid = get_grid(waves_path)

t, U_w, V_w, B_w, Bz_w, w²_w, E_w = extract_averages_timeseries(waves_directory)
t, U_s, V_s, B_s, Bz_s, w²_s, E_s = extract_averages_timeseries(stress_directory)
t, U_b, V_b, B_b, Bz_b, w²_b, E_b = extract_averages_timeseries(both_directory)

S_w = [sqrt.(U.^2 .+ V.^2) for (U, V) in zip(U_w, V_w)]
S_s = [sqrt.(U.^2 .+ V.^2) for (U, V) in zip(U_s, V_s)]
S_b = [sqrt.(U.^2 .+ V.^2) for (U, V) in zip(U_b, V_b)]

#####
##### Final profiles
#####

close("all")
fig, axs = subplots(ncols=4, sharey=true, figsize=(13, 6))

lw₁ = 2
lw₂ = 1.5
α₁ = 0.6
α₂ = 1.0

i = 1
j = 2
w_lbl = "Initial condition"

sca(axs[1])
plot(B_w[i][2:end-1], grid.zC, color=defaultcolors[j], label=w_lbl, linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[2])
plot(Bz_w[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

sca(axs[4])
plot(w²_w[i][2:end-1] / τ, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

i = length(B_w)
j = 1

sca(axs[1])
plot(B_w[i][2:end-1], grid.zC, color=defaultcolors[j], label="Growing swell, \$ t = 2 \\pi / f \$", linewidth=lw₁, alpha=α₁, linestyle="-")
plot(B_s[i][2:end-1], grid.zC, color=defaultcolors[j], label="Surface stress, no waves, \$ t = 2 \\pi / f \$", linewidth=lw₁, alpha=α₁, linestyle="--")
plot(B_b[i][2:end-1], grid.zC, color=defaultcolors[j], label="Surface stress with steady waves, \$ t = 2 \\pi / f \$", linewidth=lw₁, alpha=α₁, linestyle=":")

sca(axs[2])
plot(Bz_w[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")
plot(Bz_s[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="--")
plot(Bz_b[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle=":")

sca(axs[3])
plot(S_w[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")
plot(S_s[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="--")
plot(S_b[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle=":")

sca(axs[4])
plot(w²_w[i][2:end-1] / τ, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")
plot(w²_s[i][2:end-1] / τ, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="--")
plot(w²_b[i][2:end-1] / τ, grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle=":")

ylim(-64, 0.1)

sca(axs[1])
xlabel(L"\bar b \, \, \, (\mathrm{m \, s^{-2}})")
ylabel("\$ z \$ (m)")

#xlim(-0.0004, -0.00005)
legend(loc=3, bbox_to_anchor=(0.05, 0.02, 0.5, 0.5), prop=Dict(:size=>12), frameon=true)

sca(axs[2])

xlim(-1e-6, 4e-6)

xticks([0, 3e-6], ["0", L"3 \times 10^{-6}"])
xlabel(L"\partial_z \bar b \, \, \, (\mathrm{s^{-2}})")

sca(axs[3])
xlabel(L"\sqrt{\left ( \overline{u^\mathrm{L}} \right )^2 + \left ( \overline{v^\mathrm{L}} \right )^2} \, \, \, (\mathrm{m \, s^{-1}})")

sca(axs[4])
xlabel(L"\overline{\left (w^\mathrm{L} \right )^2} \, / \, \max(u_\star^2)")
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

savefig(joinpath(@__DIR__, "..", "figures", "figure_2.png"), dpi=480)
