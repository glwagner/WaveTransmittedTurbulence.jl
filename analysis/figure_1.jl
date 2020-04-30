using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

     run_name = "growing_waves_Qb5.0e-10_Nsq1.0e-06_f1.0e-04_dom1.5_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz384"
         part = 2
           ii = 1

run_directory = joinpath(@__DIR__, "..", "data", run_name)
     run_path = joinpath(run_directory, run_name * "_fields_part$part.jld2")
averages_path = joinpath(run_directory, run_name * "_averages.jld2")

file = jldopen(run_path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

Uˢ = uˢ(wave_amplitude, wavenumber)
τw = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
u★ = 1 #sqrt(τw) / exp(1)
grid = get_grid(run_path)

iters = get_iters(run_path)

τ(t) = τw * t / growth_time_scale * exp(-t^2 / (2 * growth_time_scale^2))

#####
##### 3D visualization
#####

t, u, v, w, b = get_fields(run_path, iters[ii])

@show f * t / 2π

umax = maximum(abs, u/u★)
@show wmax = maximum(abs, w/u★)

wlim = 0.012
wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

ulim = umax/2
ulevels = vcat([-umax], -ulim:2ulim/11:ulim, [umax])

view_elev = 50

depth = 4.0
bottom = 48

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

v_xy = v[2:end-1, 2:end-1, k+1]

close("all")

fig = figure(figsize=(14, 5))

ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4, projection="3d")
ax2 = plt.subplot2grid((1, 5), (0, 4))

plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=-1)

u★ = 1
w_im_x = ax1.contourf(w_yz / u★, YC_x, ZF_x, 
                      zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

w_im_y = ax1.contourf(XC_y, w_xz / u★, ZF_y, 
                      zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

w_im_z = ax1.contourf(XC_z, YC_z, w_xy / u★, 
                      zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

ax1.view_init(elev=view_elev, azim=-135)

cb = colorbar(w_im_x, ax=ax1, aspect=12, pad=0.05, ticks=-0.01:0.005:0.01)

α₁ = 0.6
α₂ = 0.8
lw₁ = 2
lw₂ = 1.5

ax1.set_xlabel("\$ x \$ (m)", labelpad=12.0)
ax1.set_ylabel("\$ y \$ (m)", labelpad=12.0)
ax1.set_zlabel("\$ z \$ (m)", labelpad=12.0)

ax1.set_xticks(0:48:192)
ax1.set_yticks(0:48:192)
ax1.set_xlim(0, 192)
ax1.set_ylim(0, 192)
ax1.set_zlim(-bottom, 0.1)

cb.ax.set_title(L"w \, \, \mathrm{(m \, s^{-1})}", pad=10.0)

t, U, V, B, W² = calculate_horizontal_average_timeseries(run_directory)

averages = collect_horizontal_averages(averages_path)

U = averages.U
V = averages.V
t = averages.t

@show length(t) 
i₁ = 1043
i₂ = 2086

@show f * t[i₁] / 2π
@show f * t[i₂] / 2π
@show f * (t[2]-t[1]) / 2π

bottom = 64
k_deep = searchsortedfirst(grid.zF, -bottom)
k_surface = grid.Nz

@show size(U[1]) size(grid.zC)

u★ = 1
ax2.plot(U[1, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=1, color="k",
         label=L"\langle u^\mathrm{L} \rangle \, |_{t=0}")

ax2.plot(V[1, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="--",  alpha=α₁, linewidth=1, color="k",
         label=L"\langle v^\mathrm{L} \rangle \, |_{t=0}")

ax2.plot(U[i₁, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[1], 
         label=L"\langle u^\mathrm{L} \rangle \, |_{t=\pi/2f}")

ax2.plot(V[i₁, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[1], 
         label=L"\langle v^\mathrm{L} \rangle \, |_{t=\pi/2f}")

ax2.plot(U[i₂, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="-",  alpha=α₁, linewidth=lw₁, color=defaultcolors[2], 
         label=L"\langle u^\mathrm{L} \rangle \, |_{t=\pi/f}")

ax2.plot(V[i₂, k_deep:end] / u★, grid.zC[k_deep:end], linestyle="--", alpha=α₂, linewidth=lw₂, color=defaultcolors[2], 
         label=L"\langle v^\mathrm{L} \rangle \, |_{t=\pi/f}")

legend(prop=Dict(:size=>fs), bbox_to_anchor=(0.52, 0, 0.1, 1), loc=3, markerfirst=false, frameon=false)

ax2.set_xlabel(L"\langle u^\mathrm{L} \rangle, \langle v^\mathrm{L} \rangle \, \, \mathrm{(m \, s^{-1})}", labelpad=12.0)
ax2.set_ylabel(L"z \, (\mathrm{m})", labelpad=12.0)

removespines(ax2, "right", "top")

ax3 = ax1.inset_axes([-0.05, 0.0, 0.2, 0.1])

hour = 3600
ft = range(0.0, stop=2π/f, length=1000)

ax3.plot(ft / hour, τ.(ft))
ax3.plot([1.0, 1.0] * π/f / hour, [0.0, 5e-5], "k", alpha=0.6, linewidth=1)

fs = 14
ax3.spines["right"].set_visible(false)
ax3.spines["top"].set_visible(false)

ax3.text(π/f / hour, 7.5e-5, L"\pi/f", fontsize=fs, ha="center", va="bottom")

d = 0.17
ax3.text(1.2, d, L"t", fontsize=fs, ha="center", va="center", transform=ax3.transAxes)
ax3.text(1.2, -d, "(hours)", fontsize=fs-2, ha="center", va="center", transform=ax3.transAxes)

x = -0.05
y = 1.55
δ = 0.03
ax3.text(x, y, L"\partial_t U^\mathrm{S}", fontsize=fs, ha="right", va="center", transform=ax3.transAxes)
ax3.text(x+δ, y, L"(\mathrm{m^2 \, s^{-2}})", fontsize=fs-2, ha="left", va="center", transform=ax3.transAxes)

ax3.set_xticks([0, 12])
ax3.set_xticklabels(["0", "12"], fontsize=fs)
ax3.set_yticks([0, 1e-4])
ax3.set_yticklabels(["0", L"10^{-4}"], fontsize=fs)

shift_right!(ax1, 0.02)
shift_up!(ax1, 0.01)

shift_left!(cb.ax, 0.01)
shift_up!(cb.ax, 0.12)
stretch_y!(cb.ax, -0.2)

shift_up!(ax2, 0.13)
stretch_y!(ax2, -0.2)

shift_up!(ax3, 0.06)

savefig(joinpath(@__DIR__, "..", "figures", "figure_1.png"), dpi=480)
