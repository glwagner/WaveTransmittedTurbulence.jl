using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2, FFTW, Statistics

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 16
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

function plot_w_slice!(ax, depth, path, i)

    iter = get_iters(path)[i]

    grid = get_grid(path)

    file = jldopen(path)

    f = file["coriolis/f"]
    u = file["timeseries/u/$iter"]
    v = file["timeseries/v/$iter"]
    w = file["timeseries/w/$iter"]
    t = file["timeseries/t/$iter"]

    close(file)

    @show wmax = maximum(abs, w)
    wlim = 0.5
    wlevels = vcat([-1], -wlim:2wlim/11:wlim, [1])

    #u = @. sqrt(u^2 + v^2)
    
    k = searchsortedfirst(grid.zC, -depth)

    x, y = meshgrid(grid.xC, grid.yC)

    w_xy = (
            w[2:end-1, 2:end-1, k+1] + 
            w[2:end-1, 2:end-1, k]
           ) / 2

    u_xy = (
            u[1:end-2, 2:end-1, k] + 
            u[2:end-1, 2:end-1, k]
           ) / 2

    v_xy = (
            v[2:end-1, 1:end-2, k] + 
            v[2:end-1, 2:end-1, k]
           ) / 2

    speed = @. sqrt(u_xy^2 + v_xy^2)
    max_speed = maximum(speed)

    sca(ax)
    w_im = contourf(x, y, w_xy / wmax, levels=wlevels, cmap="RdBu_r", vmin=-wlim, vmax=wlim)

    skip = 64
    hskip = skip 

    q_z = quiver(grid.Lx/2, grid.Ly/2, mean(u_xy), mean(v_xy),
                 alpha=0.7, scale=3e0 * max_speed, width=0.02)

    #q_z = quiver(   x[hskip:skip:end-hskip, hskip:skip:end-hskip],    y[hskip:skip:end-hskip, hskip:skip:end-hskip], 
    #             u_xy[hskip:skip:end-hskip, hskip:skip:end-hskip], v_xy[hskip:skip:end-hskip, hskip:skip:end-hskip],
    #             alpha=0.8, scale=5e0 * max_speed)

    text(0.05, 0.07, @sprintf("\$ t = %.2f \\times 2\\pi / f\$", f * t / 2π),
         transform=ax.transAxes, fontsize=fs,
         bbox=Dict(:alpha=>0.9, :facecolor=>"w", :edgecolor=>"w"))
    
    return w_im
end

function plot_w_spectra!(ax, depth, path, i)

    iter = get_iters(path)[i]

    grid = get_grid(path)

    file = jldopen(path)

    f = file["coriolis/f"]
    w = file["timeseries/w/$iter"]
    t = file["timeseries/t/$iter"]

    close(file)

    k = searchsortedfirst(grid.zF, -depth)
    w_xy = w[2:end-1, 2:end-1, k+1]
    
    what = fft(w_xy .+ 0im)

    what = log10.(abs.(what))

    @show wmax = maximum(what)
    @show wmin = minimum(what)

    wlim = 0 #0.2 * wmax
    wmin = -1
    wmax = 1
    wlevels = vcat([-16], wmin:(wmax-wmin)/5:wmax, [maximum(what)])

    k = fftfreq(grid.Nx, 4π/grid.Δx)
    l = fftfreq(grid.Ny, 4π/grid.Δy)
    K, L = meshgrid(k, l)

    K = fftshift(K)
    L = fftshift(L)
    what = fftshift(what)

    sca(ax)
    what_im = contourf(K, L, what, vmin=wmin, vmax=wmax, levels=wlevels)

    @show maximum(abs, what)

    return what_im
end

close("all")
fig, axs = subplots(nrows=2, ncols=6, figsize=(23, 9))
 
name = "surface_stress_with_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
directory = joinpath(@__DIR__, "..", "data", name)

paths = [
         joinpath(directory, name * "_fields_part1.jld2"),
         joinpath(directory, name * "_fields_part1.jld2"),
         joinpath(directory, name * "_fields_part2.jld2"),
         joinpath(directory, name * "_fields_part2.jld2"),
         joinpath(directory, name * "_fields_part2.jld2"),
         joinpath(directory, name * "_fields_part3.jld2"),
         #joinpath(directory, name * "_fields_part3.jld2"),
         #joinpath(directory, name * "_fields_part3.jld2"),
        ]

ii = [
      2,
      3,
      1, 
      2, 
      3, 
      1,
      #2,
      #3,
     ]

depth = 4

for (j, (i, path)) in enumerate(zip(ii, paths))
    ax = axs[1, j]
    w_im = plot_w_slice!(ax, 2, path, i)

    ax = axs[2, j]
    w_im = plot_w_slice!(ax, 8, path, i)
end

for i = 1:size(axs, 1)
    for j = 2:size(axs, 2)
        axs[i, j].tick_params(left=false, labelleft=false)
    end

    sca(axs[i, 1])
    ylabel(L"y \, \, \mathrm{(m)}")
end

for j = 1:size(axs, 2)
    axs[1, j].tick_params(bottom=false, labelbottom=false)
    sca(axs[2, j])
    xlabel(L"x \, \, \mathrm{(m)}")
end

for i = 1:size(axs, 1)
    for j = 1:size(axs, 2)
        #shift_left!(axs[i, j], 0.01)
        # stretch_x!(axs[i, j], 0.02)
        # stretch_y!(axs[i, j], 0.02)
        removespines(axs[i, j], "top", "bottom", "left", "right")
    end
end

for j = 1:size(axs, 2)
    shift_down!(axs[1, j], 0.06)
    shift_up!(axs[2, j], 0.06)
end

for i = 1:size(axs, 1), j = 1:size(axs, 2)
    axs[i, j].set_aspect(1)
end

for i = 1:size(axs, 1)
    for j = 2:size(axs, 2)
        shift_left!(axs[i, j], 0.015 * (j-1))
    end
end

δ = 0.05
for (i, depth) = zip((1, 2), (2, 8))

    depth_label = @sprintf("\$ z = - %d \$ m", depth)
    text(1.1, 0.5 + δ, L"w / \max | w |", transform=axs[i, size(axs, 2)].transAxes,
         fontsize=fs, va="center", ha="left")
    
    text(1.1, 0.5 - δ, depth_label, transform=axs[i, size(axs, 2)].transAxes,
         fontsize=fs, va="center", ha="left")
end

savefig(joinpath(@__DIR__, "..", "figures", "langmuir_rotation.png"), dpi=480)

#cb = colorbar(w_im, ax=axs, aspect=20, pad=0.10, shrink=0.8, ticks=-0.5:0.1:0.5)
#cb.ax.set_title(L"w^\mathrm{L} / \max | w^\mathrm{L} |", pad=12.0)
