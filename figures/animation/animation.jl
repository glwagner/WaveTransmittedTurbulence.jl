using WaveTransmittedTurbulence, JLD2, PyPlot, PyCall, Printf, Statistics, Glob

using Oceananigans.Buoyancy: g_Earth

second = 1.0
minute = 60second
hour = 60minute

mplot3d = pyimport("mpl_toolkits.mplot3d")

function plot_3d!(fig, grid, path, i, u★)
    @show i

    file = jldopen(path)

    w_yz = file["timeseries/w_yz/$i"]
    w_xz = file["timeseries/w_xz/$i"]
    w_xy = file["timeseries/w_xy/$i"]
       t = file["timeseries/t/$i"]

    close(file)

    wmax = max(
               maximum(abs, w_yz),
               maximum(abs, w_xz),
               maximum(abs, w_xy)
              ) / u★

    wrms = sqrt(mean([
                      mean(w_yz.^2),
                      mean(w_xz.^2),
                      mean(w_xy.^2)
                     ]) ) / u★

    @show wmax wrms
    wlim = (3 * wrms + wmax / 2) / 2
    wlim = ceil(Int, wlim * 4) / 4

    if wlim > wmax
        wlevels = vcat(-wlim:2wlim/11:wlim)
    else
        wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])
    end

    dtick = round(Int, 2 * wlim / 11 * 100) / 100
    ticks = -5dtick:dtick:5dtick

    view_elev = 50
    depth = 2.0
    bottom = 48

    k = searchsortedfirst(grid.zF, -depth)
    k_deep = searchsortedfirst(grid.zF, -bottom)

    x_offset, y_offset, z_offset = grid.Δx, grid.Δy, grid.zF[k]

    YC_x, ZF_x = meshgrid(grid.yC, grid.zF[k_deep:k])
    XC_y, ZF_y = meshgrid(grid.xC, grid.zF[k_deep:k])
    XC_z, YC_z = meshgrid(grid.xC, grid.yC)

    w_xy = w_xy[2:end-1, 2:end-1, 1]
    w_xz = w_xz[2:end-1, 1, k_deep+1:k+1]
    w_yz = w_yz[1, 2:end-1, k_deep+1:k+1]

    ax = gca()
    sca(ax)

    ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1, projection="3d")

    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, hspace=-1)

    w_im_x = ax1.contourf(w_yz / u★, YC_x, ZF_x, 
                          zdir="x", levels=wlevels, offset=x_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    w_im_y = ax1.contourf(XC_y, w_xz / u★, ZF_y, 
                          zdir="y", levels=wlevels, offset=y_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    w_im_z = ax1.contourf(XC_z, YC_z, w_xy / u★, 
                          zdir="z", levels=wlevels, offset=z_offset, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    ax1.view_init(elev=view_elev, azim=-135)

    t_hours = floor(Int, t/hour)
    t_minutes = floor(Int, (t - t_hours * hour) / minute)
    t_seconds = round(Int, t - t_hours * hour - t_minutes * minute)

    time_stamp = @sprintf("\$ t = %02d \\mathrm{h \\,} %02d \\mathrm{m \\,} %02d \\mathrm{s}  \$", t_hours, t_minutes, t_seconds)
    text2D(0.1, 0.9, time_stamp, fontsize=16, transform=ax1.transAxes)

    cb = colorbar(w_im_x, ax=ax1, ticks=ticks)

    cb.ax.set_title(L"w \, /\max(u_\star)", pad=10.0)

    ax1.set_xlabel("\$ x \$ (m)", labelpad=12.0)
    ax1.set_ylabel("\$ y \$ (m)", labelpad=12.0)
    ax1.set_zlabel("\$ z \$ (m)", labelpad=12.0)

    ax1.set_xlim(-8, 136)
    ax1.set_ylim(-8, 136)
    ax1.set_zlim(-48, 0)

    pos_cb = get_position(cb.ax)

    return nothing
end 

#####
##### A bit of setup
#####

function get_paths(case)
         run_name = "$(case)_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
    run_directory = joinpath(@__DIR__, "..", "..", "data", run_name)
           nparts = length(cd(() -> glob(run_name * "_slices*"), run_directory))
            paths = [joinpath(run_directory, run_name * "_slices_part$i.jld2") for i = 1:nparts]

    return paths
end

# Get friction velocity
growing_paths = get_paths("growing_waves")

file = jldopen(growing_paths[1])

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]

close(file)

Uˢ = uˢ(wave_amplitude, wavenumber)
 τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
u★ = sqrt(τ)

fs = 14
#plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("font"; size=fs)
#plt.rc("text", usetex=true)

#####
##### Make the animation
#####

case = "surface_stress_with_waves"
paths = get_paths(case)
           
grid = get_grid(paths[1])

close("all")
fig = figure(figsize=(14, 7))

dpi = 480

images_directory = joinpath(@__DIR__, "frames", case)
mkpath(images_directory)

@show save_points = [length(get_iters(path)) for path in paths]

# Animate
ipath = 7
path = paths[ipath]

file_iters = get_iters(path)
cumulative_images = ipath == 0 ? 0 : sum(save_points[1:ipath-1])

for (j_img, iter) in enumerate(file_iters)
    plot_3d!(fig, grid, path, iter, u★)
    image_path = joinpath(images_directory, @sprintf("%s_%06d.png", case, j_img + cumulative_images))
    savefig(image_path, dpi=dpi)
end
