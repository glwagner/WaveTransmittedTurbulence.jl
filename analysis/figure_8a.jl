include("plot_initial_condition_study_averages.jl")

axs = plot_initial_condition_study_profiles("2.0x", "Nh256_Nz256")
#savefig(joinpath(@__DIR__, "..", "figures", "figure_8.png"), dpi=480)
#
wave_amplitude = "2.0x"
suffix = "Nh256_Nz384"
resting_name = "initial_condition_study_resting_" * wave_amplitude * "_" * suffix
resting_path = joinpath(@__DIR__, "..", "data", resting_name, resting_name * "_averages.jld2")
grid = get_grid(resting_path)
resting_averages = collect_horizontal_averages(resting_path)

i = 105

sca(axs[1])
plot(resting_averages.b[i, :], grid.zC, "k")

sca(axs[2])
plot(resting_averages.bz[i, 2:end-1], grid.zF[2:end-1], "k")

sca(axs[3])
plot(resting_averages.S[i, :], grid.zC, "k")

sca(axs[4])
plot(resting_averages.W²[i, :], grid.zF, "k")
