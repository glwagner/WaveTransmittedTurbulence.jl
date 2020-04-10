include("cuda_utils.jl")
select_device!(1)

makeplot = false

include("ics_setup.jl")

Nh, Nz = has_cuda() ?
    (256, 256) :
    (64, 64)

#model = initial_conditions_model("resting", 0.8, Nh, Nz)
#run_model!(model, tf)

#model = initial_conditions_model("excited", 0.8, Nh, Nz)
#run_model!(model, tf)

model = initial_conditions_model("resting", 1.6, Nh, Nz)
run_model!(model, tf)

#model = initial_conditions_model("excited", 1.6, Nh, Nz)
#run_model!(model, tf)

#model = initial_conditions_model("control", 0.0, Nh, Nz)
#run_model!(model, wizard, tf)


exit()
