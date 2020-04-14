using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

mplot3d = pyimport("mpl_toolkits.mplot3d")

     run_name = "free_convection_Qb5.0e-10_Nsq1.0e-06_stop0.5_Nh128_Nz128"
run_directory = joinpath(@__DIR__, "..", "data", run_name)
     run_path = joinpath(run_directory, run_name * "_fields.jld2")

file = jldopen(run_path)

f = file["coriolis/f"]

close(file)

grid = get_grid(run_path)

iters = get_iters(run_path)

t, U, V, S, B, Bz, w², E = extract_averages_timeseries(run_directory, part=1)

#####
##### Final profiles
#####


fig, axs = subplots(ncols=4, sharey=true, figsize=(10, 4))

lw₁ = 2
lw₂ = 1.5
α₁ = 0.6
α₂ = 1.0

#i = length(U)
for i = 1:length(U)

    j = i

    sca(axs[1])
    plot(B[i][2:end-1], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    sca(axs[2])
    plot(Bz[i][2:end-1], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    sca(axs[3])
    plot(w²[i][2:end-1], grid.zF, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")

    sca(axs[4])
    plot(E[i], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₁, linestyle="-")
end
