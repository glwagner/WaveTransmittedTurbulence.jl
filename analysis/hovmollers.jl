using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

fs = 12
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

name = "surface_stress_with_waves_Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh128_Nz128"

directory = joinpath(@__DIR__, "..", "data", name)
     path = joinpath(directory, name * "_averages.jld2")

file = jldopen(path)

   wave_amplitude = file["surface_waves/wave_amplitude"]
       wavenumber = file["surface_waves/wavenumber"]
growth_time_scale = file["surface_waves/growth_time_scale"]
                f = file["coriolis/f"]

close(file)

averages = collect_horizontal_averages(path)
grid = get_grid(path)
iters = get_iters(path)

speed = @. sqrt(averages.U^2 + averages.V^2)
stress = @. sqrt((averages.wu + averages.τ₁₃)^2 + (averages.wv + averages.τ₂₃)^2)

TC, ZC = meshgrid(averages.t * f / 2π, grid.zC)
TF, ZF = meshgrid(averages.t * f / 2π, grid.zF)

close("all")
fig, axs = subplots(nrows=3, figsize=(24, 12), sharey=true, sharex=true)

sca(axs[1])

speed_limits = (0.03, 0.06)
levels = vcat([0.0], speed_limits[1]:0.005:speed_limits[2], [maximum(speed)])

contourf(TC, ZC, speed', vmin=speed_limits[1], vmax=speed_limits[2], levels=levels, cmap="YlGnBu_r")

sca(axs[2])

stress_limits = (1e-7, 4e-5)
levels = vcat([0.0], stress_limits[1]:1e-6:stress_limits[2], [maximum(stress)])

contourf(TF, ZF, stress', vmin=stress_limits[1], vmax=stress_limits[2], levels=levels, cmap="YlGnBu_r")

sca(axs[3])

energy_limits = (1e-6, 2e-4)
levels = vcat([0.0], energy_limits[1]:1e-5:energy_limits[2], [maximum(averages.E)])

max_cubed = maximum(abs, averages.W³)
vertical_limits = [-1, 1] * 5e-7
levels = vcat(-[max_cubed], vertical_limits[1]:1e-8:vertical_limits[2], [max_cubed])
contourf(TF, ZF, averages.W³', cmap="RdBu_r", levels=levels, vmin=vertical_limits[1], vmax=vertical_limits[2])

ylim(-48, 0)
xlim(0, 2)
