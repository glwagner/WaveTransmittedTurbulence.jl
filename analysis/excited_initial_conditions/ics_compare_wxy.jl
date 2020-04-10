include("visualization.jl")
include("analysis.jl")
include("files.jl")

using Printf

#filename = "data/initial_conditions_resting_4x_Nh192_Nz384_fields_part9.jld2"
suffix = "1x"
control_filename = "data/ics_control_Nh256_Nz256_fields_part9.jld2"
resting_filename = @sprintf("data/ics_resting_%s_Nh256_Nz256_fields_part9.jld2", suffix)
excited_filename = @sprintf("data/ics_excited_%s_Nh256_Nz256_fields_part9.jld2", suffix)

f = get_parameter(control_filename, "coriolis", "f")

grid = get_grid(resting_filename)

control_iters = get_iters(control_filename)
resting_iters = get_iters(resting_filename)
excited_iters = get_iters(excited_filename)

τ = get_wind_stress(resting_filename)

i_c = control_iters[end]
i_r = resting_iters[end]
i_e = excited_iters[end]

t_c, u_c, v_c, w_c, b_c, ν_c, κ_c = get_fields(control_filename, i_c)
t_r, u_r, v_r, w_r, b_r, ν_r, κ_r = get_fields(resting_filename, i_r)
t_e, u_e, v_e, w_e, b_e, ν_e, κ_e = get_fields(excited_filename, i_e)

b′_c, U_c, V_c, B_c, Bz_c, E_c, w²_c, e_c = calculate_statistics(u_c, v_c, w_c, b_c)
b′_r, U_r, V_r, B_r, Bz_r, E_r, w²_r, e_r = calculate_statistics(u_r, v_r, w_r, b_r)
b′_e, U_e, V_e, B_e, Bz_e, E_e, w²_e, e_e = calculate_statistics(u_e, v_e, w_e, b_e)

@show k_max_w² = argmax(w²_r)

wlim = 3
wlevels = vcat([-2wlim], -wlim:2wlim/11:wlim, [2wlim])

xCˣʸ = repeat(grid.xC, 1, grid.Ny)
yCˣʸ = repeat(reshape(grid.yC, 1, grid.Ny), grid.Nx, 1)

close("all")

fig, axs = subplots(ncols=3, figsize=(12, 4))

hpad = 0.06
fs = 12
box_props = Dict(:alpha=>0.9, :facecolor=>"w", :edgecolor=>"w")

hplot = 4
k = searchsortedfirst(grid.zF, -hplot)

sca(axs[1])
im_c = contourf(xCˣʸ, yCˣʸ, w_c[2:end-1, 2:end-1, k] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

#text(hpad*grid.Lx, hpad*grid.Lx, "Reference", fontsize=fs, bbox=box_props)
xlabel("\$ x \$ (m)")
ylabel("\$ y \$ (m)")

sca(axs[2])
im_e = contourf(xCˣʸ, yCˣʸ, w_e[2:end-1, 2:end-1, k] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

#text(hpad*grid.Lx, hpad*grid.Lx, "Excited", fontsize=fs, bbox=box_props)
xlabel("\$ x \$ (m)")

sca(axs[3])
im_r = contourf(xCˣʸ, yCˣʸ, w_r[2:end-1, 2:end-1, k] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

#text(hpad*grid.Lx, hpad*grid.Lx, "Equilibrated", fontsize=fs, bbox=box_props)
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
cbar.ax.set_title(L"w / u_\star")

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

savefig(@sprintf("w_xy_%s_h%d.png", suffix, hplot), dpi=480)
