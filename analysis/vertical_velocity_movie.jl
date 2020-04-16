using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

     run_name = "growing_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh32_Nz32"
run_directory = joinpath(@__DIR__, "..", "data", run_name)
     run_path = joinpath(run_directory, run_name * "_slices.jld2")

file = jldopen(run_path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

Uˢ = uˢ(wave_amplitude, wavenumber)
 τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
u★ = sqrt(τ)

grid = get_grid(run_path)
iters = get_iters(run_path)

function plot_3d!(fig, path, i)
    @show i

    file = jldopen(path)

    w_yz = file["timeseries/w_yz/$i"]
    w_xz = file["timeseries/w_xz/$i"]
    w_xy = file["timeseries/w_xy/$i"]

    close(file)

    wmax = max(
               maximum(abs, w_yz),
               maximum(abs, w_xz),
               maximum(abs, w_xy)
              ) / u★

    wlim = 0.8 * wmax

    wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

    view_elev = 50
    depth = 2.0
    bottom = 48

    k = searchsortedfirst(grid.zF, -depth)
    k_deep = searchsortedfirst(grid.zF, -bottom)

    x_offset, y_offset, z_offset = 0, 0, grid.zF[k]

    YC_x, ZF_x = meshgrid(grid.yC, grid.zF[k_deep:k])
    XC_y, ZF_y = meshgrid(grid.xC, grid.zF[k_deep:k])
    XC_z, YC_z = meshgrid(grid.xC, grid.yC)

    w_xy = w_xy[2:end-1, 2:end-1, 1]
    w_xz = w_xz[2:end-1, 1, k_deep+1:k+1]
    w_yz = w_yz[1, 2:end-1, k_deep+1:k+1]

    ax = gca()
    sca(ax)
    ax1 = plt.subplot2grid((1, 5), (0, 0), colspan=4, projection="3d")

    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=-1)

    w_im_x = ax1.contourf(w_yz / u★, YC_x, ZF_x, 
                          zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    w_im_y = ax1.contourf(XC_y, w_xz / u★, ZF_y, 
                          zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    w_im_z = ax1.contourf(XC_z, YC_z, w_xy / u★, 
                          zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    ax1.view_init(elev=view_elev, azim=-135)

    cb = colorbar(w_im_x, ax=ax1, ticks=wlim:0.4:wlim)

    pos_cb = get_position(cb.ax)

    return nothing
end 

close("all")
fig = figure(figsize=(14, 5))

for (j, i) in enumerate(iters)
    plot_3d!(fig, run_path, i)
    savefig("three_d_$j.png", dpi=480)
end

