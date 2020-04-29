using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2, Printf

fs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
make_axes_locatable = axes_grid1.make_axes_locatable
mplot3d = pyimport("mpl_toolkits.mplot3d")

wave_amplitude = "1.0x"
suffix = "Nh256_Nz256"

resting_name = "initial_condition_study_resting_" * wave_amplitude * "_" * suffix
excited_name = "initial_condition_study_excited_" * wave_amplitude * "_" * suffix

resting_path = joinpath(@__DIR__, "..", "data", resting_name, resting_name * "_fields_part1.jld2")
excited_path = joinpath(@__DIR__, "..", "data", excited_name, excited_name * "_fields_part1.jld2")

f = get_parameter(resting_path, "coriolis", "f")
τ = get_parameter(resting_path, "boundary_conditions", "Qᵘ")
u★ = sqrt(abs(τ))

grid = get_grid(resting_path)
iters_r = get_iters(resting_path)
iters_e = get_iters(excited_path)

t_r, u_r, v_r, w_r, b_r = get_fields(resting_path, iters_r[3])
t_e, u_e, v_e, w_e, b_e = get_fields(excited_path, iters_e[3])

@show t_r * f / 2π
@show t_e * f / 2π

wlim = 5
wlevels = vcat([-2wlim], -wlim:2wlim/11:wlim, [2wlim])

view_elev = 45
depth = 4.0
bottom = 32

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

fig = figure(figsize=(13, 5))

ax_r = plt.subplot2grid((7, 2), (1, 0), rowspan=5, projection="3d")
ax_e = plt.subplot2grid((7, 2), (1, 1), rowspan=5, projection="3d")

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
cb.ax.set_title(L"w / u_\star", pad=12.0)

ax_r.set_zlim(-bottom, 0)
ax_e.set_zlim(-bottom, 0)

ax_r.set_xlim(0, 128)
ax_r.set_ylim(0, 128)
ax_e.set_xlim(0, 128)
ax_e.set_ylim(0, 128)

ax_r.xaxis.set_ticks(0:32:128)
ax_r.yaxis.set_ticks(0:32:128)
ax_r.zaxis.set_ticks(-28:12:0)

ax_e.xaxis.set_ticks(0:32:128)
ax_e.yaxis.set_ticks(0:32:128)
ax_e.zaxis.set_ticks(-28:12:0)

ax_r.set_xlabel("\$ x \$ (m)", labelpad=24.0)
ax_r.set_ylabel("\$ y \$ (m)", labelpad=24.0)
ax_r.set_zlabel("\$ z \$ (m)", labelpad=4.0)

ax_e.set_xlabel("\$ x \$ (m)", labelpad=24.0)
ax_e.set_ylabel("\$ y \$ (m)", labelpad=24.0)
ax_e.set_zlabel("\$ z \$ (m)", labelpad=4.0)

x, y = 0.15, 0.9
ax_r.text2D(x, y, "Resting", fontsize=fs, transform=ax_r.transAxes)
ax_e.text2D(x, y, "Excited", fontsize=fs, transform=ax_e.transAxes)

shift_down!(ax_r,  0.05)
shift_down!(ax_e,  0.05)

shift_left!(ax_r, 0.03)
shift_right!(ax_e, 0.04)
shift_right!(cb.ax, 0.06)

stretch_x!(ax_r, 0.13)
stretch_x!(ax_e, 0.13)
shift_right!(cb.ax, 0.08)

stretch_y!(ax_r, 0.1)
stretch_y!(ax_e, 0.1)


savefig(joinpath(@__DIR__, "..", "figures", "figure_9.png"), dpi=480)
