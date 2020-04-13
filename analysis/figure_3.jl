using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

fs = 16
lblfs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

suffix = "Qb1.0e-09_Nsq1.0e-06_init0.3_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
part = 2

   waves_name = "growing_waves_$suffix"
  stress_name = "surface_stress_no_waves_$suffix"
    both_name = "surface_stress_with_waves_$suffix"

function name_to_path(name, part)
    directory = joinpath(@__DIR__, "..", "data", name)
    return joinpath(directory, name * "_fields_part$part.jld2"), directory
end

 waves_path,  waves_directory = name_to_path(waves_name, part)
stress_path, stress_directory = name_to_path(stress_name, part)
  both_path,   both_directory = name_to_path(both_name, part)

grid = get_grid(waves_path)
wave_amplitude, wavenumber, growth_time_scale = get_multiple_parameters(waves_path, "surface_waves", "wave_amplitude", "wavenumber", "growth_time_scale")
f = get_parameter(waves_path, "coriolis", "f")

Uˢ = uˢ(wave_amplitude, wavenumber)
 τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
u★ = sqrt(τ)

t, U_w, V_w, S_w, B_w, Bz_w, w²_w, E_w = extract_averages_timeseries(waves_directory)
t, U_s, V_s, S_s, B_s, Bz_s, w²_s, E_s = extract_averages_timeseries(stress_directory)
t, U_b, V_b, S_b, B_b, Bz_b, w²_b, E_b = extract_averages_timeseries(both_directory)

#####
##### Final profiles
#####

close("all")
fig, axs = subplots(ncols=4, sharey=true, figsize=(10, 6))

lw₁ = 2
α₁ = 0.6

plot_profiles!(axs, grid, B_w[1], Bz_w[1], nothing, w²_w[1] / τ; label="Initial condition", color=defaultcolors[2], linewidth=lw₁, alpha=α₁, linestyle="-")

i = length(B_w)
kwargs = (linewidth=lw₁, alpha=α₁, color=defaultcolors[1])
time_label = "\$ t = 2 \\pi / f \$"

plot_profiles!(axs, grid, B_w[i], Bz_w[i], S_w[i], w²_w[i] / τ; label="Growing swell, \$ t = 2 \\pi / f \$",                    linestyle="-", kwargs...)
plot_profiles!(axs, grid, B_s[i], Bz_s[i], S_s[i], w²_s[i] / τ; label="Surface stress, no waves, \$ t = 2 \\pi / f \$",         linestyle="--", kwargs...)
plot_profiles!(axs, grid, B_b[i], Bz_b[i], S_b[i], w²_b[i] / τ; label="Surface stress with steady waves, \$ t = 2 \\pi / f \$", linestyle=":", kwargs...)

ylim(-56, 0.1)

sca(axs[1])
xlabel(L"\bar b \, \, \, (\mathrm{m \, s^{-2}})", fontsize=lblfs, labelpad=16.0)
ylabel("\$ z \$ (m)", fontsize=lblfs)

xlim(-6e-5, 0.0)
xticks([-5e-5, -2e-5], [L"-5 \times 10^{-5}", L"-2 \times 10^{-5}"])
legend(loc=3, bbox_to_anchor=(0.02, 0.02, 1.0, 1.0), prop=Dict(:size=>lblfs), ncol=2, frameon=true, framealpha=1.0)

sca(axs[2])

xlim(-1e-6, 4e-6)

xticks([0, 2e-6], ["0", L"2 \times 10^{-6}"])
xlabel(L"\partial_z \bar b \, \, \, (\mathrm{s^{-2}})", fontsize=lblfs, labelpad=16.0)

sca(axs[3])
xlabel(L"\sqrt{\left ( \overline{u^\mathrm{L}} \right )^2 + \left ( \overline{v^\mathrm{L}} \right )^2} \, \, \, (\mathrm{m \, s^{-1}})", fontsize=lblfs)

sca(axs[4])
xlabel(L"\overline{\left (w^\mathrm{L} \right )^2} \, / \, \max(u_\star^2)", fontsize=lblfs, labelpad=12.0)
ylabel("\$ z \$ (m)", fontsize=lblfs)

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

positions = [get_position(ax) for ax in axs]

yshift = 0.05
for (i, pos) in enumerate(positions)
    pos[2] += yshift
    axs[i].set_position(pos)
end

savefig(joinpath(@__DIR__, "..", "figures", "figure_3.png"), dpi=480)
