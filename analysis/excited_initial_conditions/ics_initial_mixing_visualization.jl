include("visualization.jl")
include("analysis.jl")
include("files.jl")

using Oceananigans, Printf

import Oceananigans: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")
Axes3D = mplot3d.Axes3D
matplotlib = pyimport("matplotlib")

r_filename = "data/ics_resting_4x_Nh256_Nz256_fields_part2.jld2"
e_filename = "data/ics_excited_4x_Nh256_Nz256_fields_part2.jld2"

f = get_parameter(r_filename, "coriolis", "f")
τ = get_wind_stress(r_filename)
u★ = sqrt(abs(τ))

grid = get_grid(r_filename)
iters_r = get_iters(r_filename)
iters_e = get_iters(e_filename)

t_r, u_r, v_r, w_r, b_r, ν_r, κ_r = get_fields(r_filename, iters_r[1])
t_e, u_e, v_e, w_e, b_e, ν_e, κ_e = get_fields(e_filename, iters_e[1])

wlim = 5
wlevels = vcat([-2wlim], -wlim:2wlim/11:wlim, [2wlim])

view_elev = 30
depth = 1.0
bottom = 28

i = 1
j = 1
k = searchsortedfirst(grid.zF, -depth)
k_deep = searchsortedfirst(grid.zF, -bottom)

x_offset, y_offset, z_offset = 0, 0, grid.zF[k]

YC_x, ZF_x = meshgrid(grid.yC, grid.zF[k_deep:k])
XC_y, ZF_y = meshgrid(grid.xC, grid.zF[k_deep:k])
XC_z, YC_z = meshgrid(grid.xC, grid.yC)

w_xy_r = w_r[2:end-1, 2:end-1, k+1]
w_xz_r = w_r[2:end-1, j, k_deep+1:k+1]
w_yz_r = w_r[i, 2:end-1, k_deep+1:k+1]

w_xy_e = w_e[2:end-1, 2:end-1, k+1]
w_xz_e = w_e[2:end-1, j, k_deep+1:k+1]
w_yz_e = w_e[i, 2:end-1, k_deep+1:k+1]

close("all")

fig = figure(figsize=(12, 6))

ax_r = plt.subplot2grid((2, 1), (0, 0), projection="3d")
ax_e = plt.subplot2grid((2, 1), (1, 0), projection="3d")

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95) #, hspace=-1)

im_wyz_r = ax_r.contourf(w_yz_r / u★, YC_x, ZF_x, zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)
im_wxz_r = ax_r.contourf(XC_y, w_xz_r / u★, ZF_y, zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)
im_wxy_r = ax_r.contourf(XC_z, YC_z, w_xy_r / u★, zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

im_wyz_e = ax_e.contourf(w_yz_e / u★, YC_x, ZF_x, zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)
im_wxz_e = ax_e.contourf(XC_y, w_xz_e / u★, ZF_y, zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)
im_wxy_e = ax_e.contourf(XC_z, YC_z, w_xy_e / u★, zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

ax_r.view_init(elev=view_elev, azim=-135)
ax_e.view_init(elev=view_elev, azim=-135)

cb = colorbar(im_wyz_r, ax=[ax_r, ax_e], ticks=-wlim:2wlim/5:wlim)
cb.ax.set_title(L"w / u_\star")

ax_r.set_zlim(-bottom, 0)
ax_e.set_zlim(-bottom, 0)

ax_r.set_xlabel("\$ x \$ (m)", labelpad=24.0)
ax_r.set_ylabel("\$ y \$ (m)", labelpad=24.0)
ax_r.set_zlabel("\$ z \$ (m)", labelpad=4.0)

ax_e.set_xlabel("\$ x \$ (m)", labelpad=24.0)
ax_e.set_ylabel("\$ y \$ (m)", labelpad=24.0)
ax_e.set_zlabel("\$ z \$ (m)", labelpad=4.0)

fs = 12
#box_props = Dict(:alpha=>0.9, :facecolor=>"w", :edgecolor=>"w")
#ax_r.text(0, 0, 0 + 2, "Equilibrated", "x", fontsize=fs, bbox=box_props)
#ax_e.text(0, 0, 0 + 2, "Excited", "x", fontsize=fs, bbox=box_props)

x, y = 0.15, 0.9
ax_r.text2D(x, y, "Equilibrated", fontsize=fs, transform=ax_r.transAxes)
ax_e.text2D(x, y, "Excited", fontsize=fs, transform=ax_e.transAxes)

savefig("ics_initial_mixing_visualization.png", dpi=480)
