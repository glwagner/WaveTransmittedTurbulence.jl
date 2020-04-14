include("plot_growing_waves_profiles.jl")

suffix = "Qb1.0e-09_Nsq1.0e-06_init0.3_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
part = 2

fig, axs = plot_growing_wave_profiles(suffix, part)

savefig(joinpath(@__DIR__, "..", "figures", "figure_3.png"), dpi=480)
