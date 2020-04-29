include("plot_growing_waves_profiles.jl")

suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a1.0_k6.3e-02_T4.0_Nh256_Nz256"

fig, axs = plot_growing_wave_profiles(suffix, i=105, ylims=(-64, 0.1), with_stress_growing_waves=true)

savefig(joinpath(@__DIR__, "..", "figures", "figure_2.png"), dpi=480)
