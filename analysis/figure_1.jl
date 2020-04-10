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

#####
##### 3D visualization
#####

t, u, v, w, b = get_fields(run_path, iters[end])

@show f * t / 2π

umax = maximum(abs, u/u★)
@show wmax = maximum(abs, w/u★)

wlim = 0.02
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

t, U, V, B, W² = extract_averages(run_directory)

i₀ = 1
i₁ = 2
i₂ = 9

@show f * t[i₁] / 2π
@show f * t[i₂] / 2π

bottom = 36

k_deep = searchsortedfirst(grid.zF, -bottom)
k_surface = grid.Nz

@show size(U[1]) size(grid.zC)

ax2.plot(U[i₀][k_deep+1:k_surface+1] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[1], label=L"\bar u / \max(u_\star) \, |_{t=0}")
ax2.plot(V[i₀][k_deep+1:k_surface+1] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[1], label=L"\bar v / \max(u_\star) \, |_{t=0}")

ax2.plot(U[i₁][k_deep+1:k_surface+1] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[2], label=L"\bar u / \max(u_\star) \, |_{t=\pi/4f}")
ax2.plot(V[i₁][k_deep+1:k_surface+1] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[2], label=L"\bar v / \max(u_\star) \, |_{t=\pi/4f}")

ax2.plot(U[i₂][k_deep+1:k_surface+1] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[3], label=L"\bar u / \max(u_\star) \, |_{t=2\pi/f}")
ax2.plot(V[i₂][k_deep+1:k_surface+1] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[3], label=L"\bar v / \max(u_\star) \, |_{t=2\pi/f}")

legend(prop=Dict(:size=>12), bbox_to_anchor=(0.3, 0, 0.1, 1), loc=3) #, frameon=true)

ax2.set_xlabel(L"(\bar u, \bar v) / \max(u_\star)")
ax2.set_ylabel(L"z \, (\mathrm{m})", labelpad=12.0)

removespines(ax2, "right", "top")

savefig("figure_1.png", dpi=480)
