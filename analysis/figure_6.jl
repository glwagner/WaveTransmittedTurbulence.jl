using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2, Printf

fs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
make_axes_locatable = axes_grid1.make_axes_locatable

wave_amplitude = "1.0x"
suffix = "Nh256_Nz256"

control_name = "initial_condition_study_resting_0.0x_" * suffix
resting_name = "initial_condition_study_resting_" * wave_amplitude * "_" * suffix
excited_name = "initial_condition_study_excited_" * wave_amplitude * "_" * suffix

control_path = joinpath(@__DIR__, "..", "data", control_name, control_name * "_fields_part3.jld2")
resting_path = joinpath(@__DIR__, "..", "data", resting_name, resting_name * "_fields_part3.jld2")
excited_path = joinpath(@__DIR__, "..", "data", excited_name, excited_name * "_fields_part3.jld2")

f = get_parameter(control_path, "coriolis", "f")

grid = get_grid(resting_path)

control_iters = get_iters(control_path)
resting_iters = get_iters(resting_path)
excited_iters = get_iters(excited_path)

file = jldopen(resting_path)

τ = abs(file["boundary_conditions/Qᵘ"])

close(file)

i_c = control_iters[end]
i_r = resting_iters[end]
i_e = excited_iters[end]

#i_c = control_iters[2]
#i_r = resting_iters[2]
#i_e = excited_iters[2]

t_c, u_c, v_c, w_c, b_c = get_fields(control_path, i_c)
t_r, u_r, v_r, w_r, b_r = get_fields(resting_path, i_r)
t_e, u_e, v_e, w_e, b_e = get_fields(excited_path, i_e)

@show f * t_c / 2π

#b′_c, U_c, V_c, B_c, Bz_c, E_c, w²_c, e_c = calculate_statistics(u_c, v_c, w_c, b_c)
#b′_r, U_r, V_r, B_r, Bz_r, E_r, w²_r, e_r = calculate_statistics(u_r, v_r, w_r, b_r)
#b′_e, U_e, V_e, B_e, Bz_e, E_e, w²_e, e_e = calculate_statistics(u_e, v_e, w_e, b_e)

#@show k_max_w² = argmax(w²_r)

wlim = 3
wlevels = vcat([-2wlim], -wlim:2wlim/11:wlim, [2wlim])

xCˣʸ = repeat(grid.xC, 1, grid.Ny)
yCˣʸ = repeat(reshape(grid.yC, 1, grid.Ny), grid.Nx, 1)

close("all")

fig, axs = subplots(ncols=3, figsize=(12, 4))

hpad = 0.04
box_props = Dict(:alpha=>0.9, :facecolor=>"w", :edgecolor=>"w")

hplot = 4
k = searchsortedfirst(grid.zF, -hplot)

sca(axs[1])
im_c = contourf(xCˣʸ, yCˣʸ, w_c[2:end-1, 2:end-1, k] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

text(hpad*grid.Lx, hpad*grid.Lx, "Reference", ha="left", va="bottom", fontsize=fs, bbox=box_props)
xlabel("\$ x \$ (m)")
ylabel("\$ y \$ (m)")

sca(axs[2])
im_e = contourf(xCˣʸ, yCˣʸ, w_e[2:end-1, 2:end-1, k] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

text(hpad*grid.Lx, hpad*grid.Lx, "Excited", ha="left", va="bottom", fontsize=fs, bbox=box_props)
xlabel("\$ x \$ (m)")

sca(axs[3])
im_r = contourf(xCˣʸ, yCˣʸ, w_r[2:end-1, 2:end-1, k] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

text(hpad*grid.Lx, hpad*grid.Lx, "Resting", ha="left", va="bottom", fontsize=fs, bbox=box_props)
xlabel("\$ x \$ (m)")

axs[2].tick_params(left=false, labelleft=false)
axs[3].tick_params(left=false, labelleft=false)

for ax in axs
    ax.set_aspect(1, adjustable="box")
end

divider = make_axes_locatable(axs[1])
cax_c = divider.append_axes("right", size="5%", pad=0.15)

divider = make_axes_locatable(axs[2])
cax_e = divider.append_axes("right", size="5%", pad=0.15)

divider = make_axes_locatable(axs[3])
cax = divider.append_axes("right", size="5%", pad=0.15)
cbar = colorbar(im_r, cax=cax, ticks=-wlim:2wlim/5:wlim)
cbar.ax.set_title(L"w / u_\star", pad=12.0)

axs[1].yaxis.set_ticks(32:32:96)
axs[1].xaxis.set_ticks(32:32:96)
axs[2].xaxis.set_ticks(32:32:96)
axs[3].xaxis.set_ticks(32:32:96)

pos1 = [b for b in axs[1].get_position().bounds]
pos2 = [b for b in axs[2].get_position().bounds]
pos3 = [b for b in axs[3].get_position().bounds]

d = pos2[1] - pos1[1]

pos1[1] = pos2[1] - 0.9*d
pos3[1] = pos2[1] + 0.9*d

pos1[3] = 1.1 * pos1[3]
pos2[3] = 1.1 * pos2[3]
pos3[3] = 1.1 * pos3[3]

axs[1].set_position(pos1)
axs[2].set_position(pos2)
axs[3].set_position(pos3)

for cax in (cax_c, cax_e)
    cax.set_facecolor("None")
    for side in ("top", "bottom", "left", "right")
        cax.spines[side].set_visible(false)
    end
    cax.tick_params(left=false, labelleft=false, bottom=false, labelbottom=false)
end

savefig(joinpath(@__DIR__, "..", "figures", "figure_6.png"), dpi=480)
