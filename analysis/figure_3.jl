using WaveTransmittedTurbulence, Oceananigans, PyPlot, PyCall, Printf, JLD2

using Oceananigans.Buoyancy: g_Earth

fs = 16
lblfs = 14
plt.rc("font"; family="serif", serif=["Computer Modern Roman"], size=fs)
plt.rc("text", usetex=true)

suffix = "Qb1.0e-09_Nsq1.0e-06_init0.3_a2.0_k6.3e-02_T4.0_Nh256_Nz256"
part = 2

 waves_name = "growing_waves_$suffix"
stress_name = "surface_stress_no_waves_$suffix"
  both_name = "surface_stress_with_waves_$suffix"

include("growing_waves_profiles.jl")

savefig(joinpath(@__DIR__, "..", "figures", "figure_3.png"), dpi=480)
