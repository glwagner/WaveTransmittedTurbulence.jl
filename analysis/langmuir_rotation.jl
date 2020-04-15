using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 12
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

run_name = "surface_stress_with_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh128_Nz128"

run_directory = joinpath(@__DIR__, "..", "data", run_name)
     run_path = joinpath(run_directory, run_name * "_fields.jld2")

file = jldopen(run_path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

iters = get_iters(run_path)      

ncols = 6
nrows = ceil(Int, length(iters)/ncols)

close("all")
fig, axs = subplots(ncols=ncols, nrows=nrows)

for (iplot, iter) in enumerate(iters)

    Uˢ = uˢ(wave_amplitude, wavenumber)
     τ = wave_amplitude^2 * √(g_Earth * wavenumber) / (2 * growth_time_scale)
    u★ = sqrt(τ)
    grid = get_grid(run_path)

    t, u, v, w, b = get_fields(run_path, iter)

    @show f * t / 2π

    umax = maximum(abs, u/u★)
    @show wmax = maximum(abs, w/u★)

    wlim = 0.1 #0.5 * wmax #1 #0.3 #19 * wmax / 20#1.5
    wlevels = vcat([-wmax], -wlim:2wlim/11:wlim, [wmax])

    ulim = umax/2
    ulevels = vcat([-umax], -ulim:2ulim/11:ulim, [umax])

    depth = 2.0

    k = searchsortedfirst(grid.zF, -depth)

    x, y = meshgrid(grid.xC, grid.yC)

    w_xy = w[2:end-1, 2:end-1, k+1]
    u_xy = 0.5 * (u[2:end-1, 2:end-1, k+2] .+ u[3:end, 2:end-1, k+2])
    v_xy = 0.5 * (v[2:end-1, 2:end-1, k+2] .+ v[2:end-1, 3:end, k+2])

    @show row = ceil(Int, iplot/ncols)
    @show col = mod1(iplot, ncols)
    sca(axs[row, col])

    w_im = contourf(x, y, w_xy / u★, levels=wlevels, cmap="RdBu_r", vmin=-wlim, vmax=wlim)
                          
    skip = 20
    q_z = quiver(   x[1:skip:end, 1:skip:end],    y[1:skip:end, 1:skip:end], 
                 u_xy[1:skip:end, 1:skip:end], v_xy[1:skip:end, 1:skip:end])


    title(@sprintf("\$ t = %.2f \$", f * t / 2π))
end

for ax in axs
    ax.set_aspect(1)
end

cb = colorbar(w_im, ax=axs, ticks=-wlim:0.5:wlim)
