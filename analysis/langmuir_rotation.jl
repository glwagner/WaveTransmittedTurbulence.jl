using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 12
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

function plot_w_u!(ax, path, iter)
    grid = get_grid(path)

    file = jldopen(path)

    f = file["coriolis/f"]
    u = file["timeseries/u/$iter"]
    v = file["timeseries/v/$iter"]
    w = file["timeseries/w/$iter"]
    t = file["timeseries/t/$iter"]

    close(file)

    @show f * t / 2π
    @show wmax = maximum(abs, w)

    wlim = 0.5 * wmax
    wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

    depth = 4.0
    k = searchsortedfirst(grid.zF, -depth)

    x, y = meshgrid(grid.xC, grid.yC)

    w_xy = w[2:end-1, 2:end-1, k+1]
    u_xy = 0.5 * (u[2:end-1, 2:end-1, k+2] .+ u[3:end, 2:end-1, k+2])
    v_xy = 0.5 * (v[2:end-1, 2:end-1, k+2] .+ v[2:end-1, 3:end, k+2])

    sca(ax)
    w_im = contourf(x, y, w_xy, levels=wlevels, cmap="RdBu_r", vmin=-wlim, vmax=wlim)
                          
    skip = 20
    q_z = quiver(   x[1:skip:end, 1:skip:end],    y[1:skip:end, 1:skip:end], 
                 u_xy[1:skip:end, 1:skip:end], v_xy[1:skip:end, 1:skip:end])


    title(@sprintf("\$ t = %.2f \$", f * t / 2π))
end

name = "surface_stress_with_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
directory = joinpath(@__DIR__, "..", "data", name)

paths = [
         joinpath(directory, name * "_fields_part1.jld2"),
         joinpath(directory, name * "_fields_part2.jld2"),
         joinpath(directory, name * "_fields_part3.jld2"),
        ]

ncols = 3
nrows = length(paths)

close("all")
fig, axs = subplots(ncols=ncols, nrows=nrows, sharex=true, sharey=true)

for (path, i) = zip(paths, 1:nrows)
    iters = get_iters(path)
    for (iter, j) = zip(iters, 1:ncols)
        ax = axs[i, j]
        plot_w_u!(ax, path, iter)
    end
end

for ax in axs
    ax.set_aspect(1)
end
