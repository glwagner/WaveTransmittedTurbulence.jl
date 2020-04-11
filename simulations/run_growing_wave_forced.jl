using WaveTransmittedTurbulence

using 
    Oceananigans, 
    Oceananigans.Diagnostics, 
    Oceananigans.OutputWriters, 
    Oceananigans.Coriolis, 
    Oceananigans.BoundaryConditions, 
    Oceananigans.Forcing, 
    Oceananigans.Utils, 
    Oceananigans.TurbulenceClosures, 
    Oceananigans.Buoyancy

using Random, Printf, Glob, ArgParse

"Returns a dictionary of command line arguments."
function parse_command_line_arguments()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "--spinup"
            help = "The name of the directory to look for spinup data."
            default = "free_convection_Qb1.0e-09_Nsq1.0e-05_Nh32_Nz32"
            arg_type = String

        "--case"
            help = """The case to run. The options are:
                   1. "growing_waves"
                   2. "surface_stress_no_waves"
                   3. "surface_stress_with_waves"
                   """
            default = "growing_waves"
            arg_type = String

        "--device", "-d"
            help = "The CUDA device index on which to run the simulation."
            default = 0
            arg_type = Int
    end

    return parse_args(settings)
end

args = parse_command_line_arguments()

@hascuda select_device!(args["device"])
 
# # Set up simulation from spinup state
spinup_name = args["spinup"]

spinup_directory = joinpath(@__DIR__, "..", "data", spinup_name)

# Pick last spinup file:
filenames = cd(() -> glob("*fields*"), spinup_directory)
sortby(filename) = parse(Int, filename[length(spinup_name)+13:end-5])
sort!(filenames, by=sortby, rev=true)

filepath = joinpath(@__DIR__, "..", "data", spinup_name, filenames[1])

# Grid
grid = get_grid(filepath)

# # Stokes drift parameters
      #wave_number = 2π / 100
      wave_number = 2π / 60
   wave_amplitude = 2.0
growth_time_scale = 4hour

case = args["case"]

if case === "growing_waves"

    stokes_drift = GrowingStokesDrift(wave_number=wave_number, wave_amplitude=wave_amplitude,
                                      growth_time_scale=growth_time_scale)

    u_bcs = UVelocityBoundaryConditions(grid) # default

elseif case == "surface_stress_no_waves"

    stokes_drift = nothing

    Qᵘ = EffectiveStressGrowingStokesDrift(wave_number=wave_number, wave_amplitude=wave_amplitude,
                                           growth_time_scale=growth_time_scale)

    u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ)) # default

elseif case == "surface_stress_with_waves"

    stokes_drift = SteadyStokesDrift(wave_number=wave_number, wave_amplitude=wave_amplitude,
                                     growth_time_scale=growth_time_scale)

    Qᵘ = EffectiveStressGrowingStokesDrift(wave_number=wave_number, wave_amplitude=wave_amplitude,
                                           growth_time_scale=growth_time_scale)

    u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ)) # default

end
   
# Boundary conditions
N² = get_parameter(filepath, "initial_conditions", "N²")

b_bcs = TracerBoundaryConditions(grid, bottom=BoundaryCondition(Gradient, N²))

# Sponge layer
δ = get_parameter(filepath, "sponge_layer", "δ")
τ = get_parameter(filepath, "sponge_layer", "τ")

u_forcing = ParameterizedForcing(Fu, (δ=δ, τ=τ))
v_forcing = ParameterizedForcing(Fv, (δ=δ, τ=τ))
w_forcing = ParameterizedForcing(Fw, (δ=δ, τ=τ))
b_forcing = ParameterizedForcing(Fb, (δ=δ, τ=τ, dbdz=N²))

# Reconstruct eddy diffusivity model
Cᴬᴹᴰ = SurfaceEnhancedModelConstant(filepath)

f = get_parameter(filepath, "coriolis", "f")

# # Model instantiation and initial condition

model = IncompressibleModel(       architecture = has_cuda() ? GPU() : CPU(),
                                           grid = grid,
                                        tracers = :b,
                                       buoyancy = BuoyancyTracer(),
                                       coriolis = FPlane(f=f),
                                  surface_waves = stokes_drift,
                                        closure = AnisotropicMinimumDissipation(C=Cᴬᴹᴰ),
                            boundary_conditions = (b=b_bcs, u=u_bcs),
                                        forcing = ModelForcing(u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)
                           )

set_from_file!(model, filepath)

# # Prepare the simulation

stop_time = 4π / f

# Adaptive time-stepping
wizard = TimeStepWizard(       cfl = 0.2,
                                Δt = 1e-1,
                        max_change = 1.1,
                            max_Δt = 10.0)

messenger = SimulationProgressMessenger(model, wizard)

simulation = Simulation(model, Δt=wizard, stop_time=stop_time, progress_frequency=100, progress=messenger)

# # Specify output

prefix = @sprintf("growing_wave_forced_Qb%.1e_a%.1f_k%.1e_T%.1f_Nh%d_Nz%d",
                  get_parameter(filepath, "boundary_conditions", "Qᵇ"),
                  stokes_drift.∂z_uˢ.wave_amplitude,
                  stokes_drift.∂z_uˢ.wave_number,
                  stokes_drift.∂z_uˢ.growth_time_scale / hour,
                  model.grid.Nx, model.grid.Nz)

data_directory = joinpath(@__DIR__, "..", "data", prefix) # save data in /data/prefix

"Save a few things that we might want when we analyze the data."
function init(file, model; kwargs...)
    file["sponge_layer/δ"] = δ
    file["sponge_layer/τ"] = τ
    file["initial_conditions/N²"] = N²
    file["surface_waves/wave_number"] = model.surface_waves.∂z_uˢ.wave_number
    file["surface_waves/wave_amplitude"] = model.surface_waves.∂z_uˢ.wave_amplitude
    file["surface_waves/growth_time_scale"] = model.surface_waves.∂z_uˢ.growth_time_scale
    return nothing
end

# Three-dimensional field output
fields_to_output = merge(model.velocities, model.tracers, (νₑ=model.diffusivities.νₑ,),
                         prefix_tuple_names(:κₑ, model.diffusivities.κₑ))

field_writer = JLD2OutputWriter(model, FieldOutputs(fields_to_output); force=true, init=init,
                                    interval = π / 2f, # every quarter period
                                max_filesize = 2GiB,
                                         dir = data_directory,
                                      prefix = prefix * "_fields")

simulation.output_writers[:fields] = field_writer

#=
# Horizontal averages
averages_writer = JLD2OutputWriter(model, horizontal_averages(model); force=true, init=init,
                                   interval = 1hour, 
                                        dir = data_directory,
                                     prefix = prefix * "_averages")

simulation.output_writers[:averages] = averages_writer
=#

# # Run

print_banner(simulation)

run!(simulation)

exit() # Release GPU memory
