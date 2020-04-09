include("visualization.jl")
include("analysis.jl")
include("files.jl")

using Oceananigans, Printf

import Oceananigans: g_Earth

png_filepath = "growth_intro_plot.png"
cmap = "YlGnBu"

dir = "data"
prefix = "growth_wave_forced_convecting_spinup_ramp_a1.5_Nh256_Nz192"
fields_filename = joinpath(dir, prefix * "_fields_part5.jld2")

averages_filename = joinpath(dir, prefix * "_averages.jld2")

aˢʷ, kˢʷ = get_surface_wave_parameters(fields_filename)
T = 4hour
Qᵘ₀ = get_parameter(fields_filename, "boundary_conditions", "Qᵘ₀")
  f = get_parameter(fields_filename, "coriolis", "f")

if Qᵘ₀ === nothing
    Qᵘ₀ = -5e-6
end

 Uˢ = aˢʷ^2 * kˢʷ * sqrt(g_Earth * kˢʷ)
Qᵘ₁ = -aˢʷ^2 * sqrt(g_Earth * kˢʷ) / 2T
  τ = abs(Qᵘ₁)

  u★ = sqrt(τ)

grid = get_grid(fields_filename)

iters = get_iters(fields_filename)
last_iter = iters[end]
first_iter = iters[1]

#####
##### 3D visualization
#####

#t, u, v, w, b, ν, κ = get_fields(fields_filename, last_iter)
t, u, v, w, b, ν, κ = get_fields(fields_filename, first_iter)

@show f * t / 2π

t, U, V, B, Bz, u², v², w² = get_averages(averages_filename)

u★ = sqrt(abs(τ))

umax = maximum(abs, u/u★)
wmax = maximum(abs, w/u★)

wlim = 0.4
wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

ulim = umax/2
ulevels = vcat([-umax], -ulim:2ulim/11:ulim, [umax])

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

close("all")

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

w_im_x = ax1.contourf(w_yz / u★, YC_x, ZF_x, 
                      zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

w_im_y = ax1.contourf(XC_y, w_xz / u★, ZF_y, 
                      zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

w_im_z = ax1.contourf(XC_z, YC_z, w_xy / u★, 
                      zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

ax1.view_init(elev=view_elev, azim=-135)

cb = colorbar(w_im_x, ax=ax1, ticks=-wlim:0.1:wlim)

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

cb.ax.set_title(L"w / \max(u_\star)", pad=10.0)

i₀ = 1
i₁ = 26
i₂ = 201

@show f * t[i₁] / 2π
@show f * t[i₂] / 2π

bottom = 36

k_deep = searchsortedfirst(grid.zF, -bottom)

ax2.plot(U[i₀, k_deep:end]   / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[1], label=L"\bar u / \max(u_\star) \, |_{t=0}")
ax2.plot(V[i₀, k_deep:end]   / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[1], label=L"\bar v / \max(u_\star) \, |_{t=0}")

ax2.plot(U[i₁, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[2], label=L"\bar u / \max(u_\star) \, |_{t=\pi/4f}")
ax2.plot(V[i₁, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[2], label=L"\bar v / \max(u_\star) \, |_{t=\pi/4f}")

ax2.plot(U[i₂, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[3], label=L"\bar u / \max(u_\star) \, |_{t=2\pi/f}")
ax2.plot(V[i₂, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[3], label=L"\bar v / \max(u_\star) \, |_{t=2\pi/f}")

legend(prop=Dict(:size=>12), bbox_to_anchor=(0.3, 0, 0.1, 1), loc=3) #, frameon=true)

ax2.set_xlabel(L"(\bar u, \bar v) / \max(u_\star)")
ax2.set_ylabel(L"z \, (\mathrm{m})", labelpad=12.0)

removespines(ax2, "right", "top")
#ax2.tick_params(left=false, labelleft=false, right=true, labelright=true)
#ax2.tick_params(left=false, labelleft=false, right=true, labelright=true)
#ax2.yaxis.set_label_position("right")

savefig("growth_intro_plot.png", dpi=480)
