using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2, Printf

ticker = pyimport("matplotlib.ticker")
ScalarFormatter = ticker.ScalarFormatter

fs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

minute = 60.0
hour = 60minute

wave_amplitude = "1.0x"
suffix = "Nh256_Nz256"

control_name = "initial_condition_study_resting_0.0x_" * suffix
excited_name = "initial_condition_study_excited_" * wave_amplitude * "_" * suffix

control_path = joinpath(@__DIR__, "..", "data", control_name, control_name * "_averages.jld2")
excited_path = joinpath(@__DIR__, "..", "data", excited_name, excited_name * "_averages.jld2")

grid = get_grid(excited_path)

control_averages = collect_horizontal_averages(control_path)
excited_averages = collect_horizontal_averages(excited_path)

#####
##### Plot
#####

close("all")
fig, axs = subplots(ncols=4, figsize=(10, 4), sharey=true)

α₁ = 0.6
α₂ = 0.8

lw₂ = 2.7
lw₁ = 1.5

depth = 36

i = 210

k = searchsortedfirst(grid.zF, -depth)

@show length(control_averages.t)
@show length(excited_averages.t)

c_labels, r_labels, e_labels = Tuple([] for i = 1:3)

f = get_parameter(excited_path, "coriolis", "f")

k = 0.105
a = 0.8
g = 9.81

zC = grid.zC
@show k_c = argmax(control_averages.bz[i, :])
@show k_e = argmax(excited_averages.bz[i, :])

z_c = abs(grid.zF[k_c])
z_e = abs(grid.zF[k_e])

uS = @. a^2 * k * √(g*k) * exp(2k * zC)

c_color = "k"
e_color = "k"

sca(axs[1])
plot(control_averages.U[i, :],       grid.zC / z_c, label="Reference", color=c_color, linewidth=lw₂, alpha=α₁, "-")
plot(excited_averages.U[i, :] .- uS, grid.zC / z_e, label="Excited", color=c_color, linewidth=lw₁, alpha=α₂, "--")

sca(axs[2])
plot(control_averages.V[i, :], grid.zC / z_c, color=c_color, linewidth=lw₂, alpha=α₁, "-")
plot(excited_averages.V[i, :], grid.zC / z_e, color=c_color, linewidth=lw₁, alpha=α₂, "--")

sca(axs[3])
plot(control_averages.W²[i, :], grid.zF / z_c, color=c_color, linewidth=lw₂, alpha=α₁, "-")
plot(excited_averages.W²[i, :], grid.zF / z_e, color=c_color, linewidth=lw₁, alpha=α₂, "--")

sca(axs[4])
plot(control_averages.U[i, :], grid.zC / z_c, color=c_color, linewidth=lw₂, alpha=α₁, "-")
plot(excited_averages.U[i, :], grid.zC / z_e, color=c_color, linewidth=lw₁, alpha=α₂, "--")

ylim(-1.5, 0)

sca(axs[1])
xlabel(L"\langle u^\mathrm{E} \rangle = \langle u^\mathrm{L} \rangle - u^\mathrm{S} \, \, (\mathrm{m \, s^{-1}})", labelpad=8.0)
ylabel("\$ z / h \$")
legend(frameon=false, loc="lower right", markerfirst=false, bbox_to_anchor=(0.4, 0.01, 1.0, 1.0))

sca(axs[2])
xlabel(L"\langle v^\mathrm{L} \rangle = \langle v^\mathrm{E} \rangle \, \, (\mathrm{m \, s^{-1}})", labelpad=8.0)
legend(loc=3, prop=Dict("size" => 12), bbox_to_anchor=(0.3, 0., 1.0, 0.5), frameon=false)

sca(axs[3])
xlabel(L"\left \langle \left ( w^{\mathrm{L}} \right )^2 \right \rangle \, \, (\mathrm{m^2 \, s^{-2}})")

sca(axs[4])
xlabel(L"\langle u^\mathrm{L} \rangle \, \, (\mathrm{m \, s^{-1}})", labelpad=8.0)
ylabel("\$ z / h\$")

for (ax, lbl) in zip(axs, ("(\\textit{a})", "(\\textit{b})", "(\\textit{c})", "(\\textit{d})"))
    sca(ax)
    text(0.05, 1.03, lbl, transform=ax.transAxes, ha="left", va="bottom")
end

removespines(axs[1], "top", "right")
removespines(axs[end], "top", "left")

for j = 2:3
    removespines(axs[j], "top", "right", "left")
    axs[j].tick_params(left=false, labelleft=false)
end

axs[end].yaxis.set_label_position("right")
axs[end].tick_params(left=false, labelleft=false, right=true, labelright=true)
axs[end].set_zorder(-1)

for (i, ax) in enumerate(axs)
    stretch_x!(ax, -0.01)
    stretch_y!(ax, -0.05)
    shift_up!(ax, 0.07)
    shift_left!(ax, 0.03)
    shift_right!(ax, 0.015 * (i-1))
end

savefig(joinpath(@__DIR__, "..", "figures", "figure_10.png"), dpi=480)
