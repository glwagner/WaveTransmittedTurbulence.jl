include("plot_initial_condition_study_averages.jl")

axs = plot_initial_condition_study_profiles("1.0x", "Nh256_Nz256")
savefig(joinpath(@__DIR__, "..", "figures", "figure_7.png"), dpi=480)
