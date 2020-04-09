include("visualization.jl")
include("analysis.jl")
include("files.jl")

using Oceananigans, Printf

import Oceananigans: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")
Axes3D = mplot3d.Axes3D
matplotlib = pyimport("matplotlib")

png_filepath = "growth_intro_plot.png"
cmap = "YlGnBu"

dir = "data"
#prefix = "growth_wave_forced_strong_spinup_ramp_a1.5_Nh256_Nz192"
prefix = "spinup_convecting_only_Nh256_Nz192"
fields_filename = joinpath(dir, prefix * "_fields_part17.jld2")

averages_filename = joinpath(dir, prefix * "_averages.jld2")

grid = get_grid(fields_filename)

iters = get_iters(fields_filename)

t, u, v, w, b, ν, κ = get_fields(fields_filename, iters[1])
t, U, V, B, Bz, u², v², w² = get_averages(averages_filename)

wmax = maximum(abs, w)

@show wlim = wmax / 2
wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

view_elev = 50
depth = 2.0
bottom = 24

i = 1
j = 1
k = searchsortedfirst(grid.zF, -depth)
k_deep = searchsortedfirst(grid.zF, -bottom)

x_offset, y_offset, z_offset = 0, 0, grid.zF[k]

YC_x, ZF_x = meshgrid(grid.yC, grid.zF[k_deep:k])
XC_y, ZF_y = meshgrid(grid.xC, grid.zF[k_deep:k])
XC_z, YC_z = meshgrid(grid.xC, grid.yC)

w_xy = w[2:end-1, 2:end-1, k+1]
w_xz = w[2:end-1, j, k_deep+1:k+1]
w_yz = w[i, 2:end-1, k_deep+1:k+1]

u_xy = u[2:end-1, 2:end-1, k+1]
u_xz = u[2:end-1, j, k_deep+1:k+1]
u_yz = u[i, 2:end-1, k_deep+1:k+1]

#close("all")

fig = figure(figsize=(14, 5))

ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4, projection="3d")
ax2 = plt.subplot2grid((1, 5), (0, 4))

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=-1)

pos2 = get_position(ax2)

xshift = -0.05
yshift = 0.1
pos2[1] += xshift
pos2[2] += yshift
pos2[4] -= 2yshift

ax2.set_position(pos2)

w_im_x = ax1.contourf(w_yz, YC_x, ZF_x, 
                      zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

w_im_y = ax1.contourf(XC_y, w_xz, ZF_y, 
                      zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

w_im_z = ax1.contourf(XC_z, YC_z, w_xy, 
                      zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

ax1.view_init(elev=view_elev, azim=-135)

cb = colorbar(w_im_x, ax=ax1)

pos_cb = get_position(cb.ax)

α₁ = 0.6
α₂ = 0.8
lw₁ = 2
lw₂ = 1.5

xshift = 0.04
yshrink = 0.12

pos_cb[1] -= xshift
pos_cb[2] += yshrink
pos_cb[4] -= 2yshrink

cb.ax.set_position(pos_cb)

ax1.set_xlabel("\$ x \$ (m)", labelpad=12.0)
ax1.set_ylabel("\$ y \$ (m)", labelpad=12.0)
ax1.set_zlabel("\$ z \$ (m)", labelpad=12.0)

cb.ax.set_title(L"w", pad=10.0)

ii = [1, 11, 21, 51, 101, 401, 801]
k_deep = 1 #searchsortedfirst(grid.zF, -grid.zF[end])

for (j, i) in enumerate(ii)
    ax2.plot(Bz[i, k_deep:end], grid.zF[k_deep+1:end-1], linestyle="-", alpha=α₁, linewidth=lw₁, color=defaultcolors[j])
end

ax2.set_xlabel(L"\partial_z \bar b")
ax2.set_ylabel(L"z \, (\mathrm{m})", labelpad=12.0)
ax2.set_ylim(-16, 0)

removespines(ax2, "left", "top")
ax2.tick_params(left=false, labelleft=false, right=true, labelright=true)
ax2.yaxis.set_label_position("right")

fig, axs = subplots(ncols=2)

for (j, i) in enumerate(ii)
    sca(axs[1])
    plot(U[i, :], grid.zC, linestyle="-", alpha=α₁, linewidth=lw₁, color=defaultcolors[j])
    ylim(-16, 0)

    sca(axs[2])
    plot(V[i, :], grid.zC, linestyle="-", alpha=α₁, linewidth=lw₁, color=defaultcolors[j])
    ylim(-16, 0)
end
