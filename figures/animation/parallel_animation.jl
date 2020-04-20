@everywhere begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, "..", ".."))
    Pkg.instantiate()
end

@everywhere begin
    using WaveTransmittedTurbulence, JLD2, PyPlot, PyCall, Printf, Statistics

    using Oceananigans.Buoyancy: g_Earth

    second = 1.0
    minute = 60second
    hour = 60minute

    mplot3d = pyimport("mpl_toolkits.mplot3d")
end

@everywhere function plot_3d!(fig, grid, path, i)
    @show i

    file = jldopen(path)

    w_yz = file["timeseries/w_yz/$i"]
    w_xz = file["timeseries/w_xz/$i"]
    w_xy = file["timeseries/w_xy/$i"]
       t = file["timeseries/t/$i"]

       wave_amplitude = file["surface_waves/wave_amplitude"]
           wavenumber = file["surface_waves/wavenumber"]
    growth_time_scale = file["surface_waves/growth_time_scale"]
                    f = file["coriolis/f"]

    close(file)

    Uˢ = uˢ(wave_amplitude, wavenumber)
     τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
    u★ = 1 #sqrt(τ)

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
    wlim = 1/2 * (6 * wrms + wmax / 2)
    wlim = ceil(Int, wlim * 2) / 2

    if wlim > wmax
        wlevels = vcat(-wlim:2wlim/11:wlim)

        dtick = round(Int, 2 * wlim / 11 * 100) / 100
        ticks = -5dtick:dtick:5dtick
    else
        wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

        dtick = round(Int, 2 * wlim / 11 * 100) / 100
        ticks = vcat([-wmax], -5dtick:dtick:5dtick, [wmax])
    end

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

    pos_cb = get_position(cb.ax)

    return nothing
end 

@everywhere begin
    #####
    ##### Make the plot
    #####
    
    fs = 14
    plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
    plt.rc("text", usetex=true)

    #for case in ("growing_waves", "surface_stress_no_waves", "surface_stress_with_waves")
    case = "surface_stress_no_waves"
         run_name = "$(case)_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
    run_directory = joinpath(@__DIR__, "..", "..", "data", run_name)
            paths = [joinpath(run_directory, run_name * "_slices_part$i.jld2") for i = 1:7]
               
    grid = get_grid(paths[1])
    total_iters = 0

    close("all")
    fig = figure(figsize=(14, 7))

    dpi = 480

    images_directory = joinpath(@__DIR__, "frames", case)
    mkpath(images_directory)

    iters = [get_iters(path) for path in paths]
end

# Animate
for (ip, p) in enumerate(procs()[1:7])

    @spawnat p begin
        path = paths[ip]

        @show cumulative_iters = sum(iters[1:ip-1])
        these_iters = get_iters(path)

        for (j_img, iter) in enumerate(these_iters)
            plot_3d!(fig, grid, path, iter)
            image_path = joinpath(images_directory, @sprintf("%s_%06d.png", case, j_img + cumulative_iters))
            savefig(image_path, dpi=dpi)
        end
    end
end
