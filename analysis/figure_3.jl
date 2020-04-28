include("plot_growing_waves_profiles.jl")

#suffix = "Qb5.0e-10_Nsq1.0e-06_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
suffix = "Qb5.0e-10_Nsq1.0e-06_f1.0e-04_dom1.5_init0.5_a2.0_k6.3e-02_T4.0_Nh256_Nz384"

fig, axs = plot_growing_wave_profiles(suffix, ylims=(-64, 0.1))

savefig(joinpath(@__DIR__, "..", "figures", "figure_3.png"), dpi=480)
