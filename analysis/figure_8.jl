include("plot_initial_condition_study_averages.jl")

axs = plot_initial_condition_study_profiles("2.0x", "Nh256_Nz256")
savefig(joinpath(@__DIR__, "..", "figures", "figure_8.png"), dpi=480)
