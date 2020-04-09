include("cuda_utils.jl")
select_device!(3)

makeplot = false

include("growth_simulation_setup.jl")

#model = wave_forced_model("growth_wave_forced_convecting_spinup")
#run_model!(model, duration)

#model = wind_forced_model("growth_wind_forced_convecting_spinup")
#run_model!(model, duration)

model = wind_with_waves_model("growth_wind_with_waves_convecting_spinup")
run_model!(model, duration)

exit()
