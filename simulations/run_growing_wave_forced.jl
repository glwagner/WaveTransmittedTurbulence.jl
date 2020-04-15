#= run_growing_wave_forced.jl

This script sets up and runs a simulation of a turbulent boundary layer
forced by a growing surface wave field.

The script is intended to be used within the `WaveTransmittedTurbulence` 
environment.
=#
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
            default = "free_convection_Qb5.0e-10_Nsq1.0e-06_stop0.5_Nh256_Nz256"
            arg_type = String

        "--spinup-part"
            help = """The spinup file to look in for an initial condition.
                      0 chooses the largest part number in the spinup directory."""
            default = 0
            arg_type = Int

        "--spinup-save-point"
            help = """The save point within the indicated spinup file to use for
                      an initial condition. 0 chooses the last iteration saved in the
                      indicated spinup file."""
            default = 0
            arg_type = Int

        "--case"
            help = """The case to run. The options are:
                   1. "growing_waves"
                   2. "surface_stress_no_waves"
                   3. "surface_stress_with_waves"
                   """
            default = "growing_waves"
            arg_type = String

        "--wave-amplitude"
            help = "The equilibrium wave field amplitude in meters."
            default = 1.0
            arg_type = Float64

        "--growth-time-scale"
            help = "The time-scale for wave growth in hours."
            default = 4.0
            arg_type = Float64

        "--wavelength"
            help = "The wave length of the surface waves in meters."
            default = 100.0
            arg_type = Float64

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

# Sort spinup files so that one may be picked:
filenames = cd(() -> glob("*fields*"), spinup_directory)
sortby(filename) = parse(Int, filename[length(spinup_name)+13:end-5])
sort!(filenames, by=sortby)
part = args["spinup-part"] == 0 ? length(filenames) : args["spinup-part"]

filepath = joinpath(@__DIR__, "..", "data", spinup_name, filenames[part])

# Grid
grid = get_grid(filepath)

# # Stokes drift parameters
       wavenumber = 2π / args["wavelength"]
   wave_amplitude = args["wave-amplitude"]
growth_time_scale = args["growth-time-scale"] * hour

case = args["case"]

if case === "growing_waves"

    stokes_drift = GrowingStokesDrift(wavenumber=wavenumber, wave_amplitude=wave_amplitude,
                                      growth_time_scale=growth_time_scale)

    u_bcs = UVelocityBoundaryConditions(grid) # default

elseif case == "surface_stress_no_waves"

    stokes_drift = nothing

    Qᵘ = EffectiveStressGrowingStokesDrift(wavenumber=wavenumber, wave_amplitude=wave_amplitude,
                                           growth_time_scale=growth_time_scale)

    u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ))

elseif case == "surface_stress_with_waves"

    stokes_drift = SteadyStokesDrift(wavenumber=wavenumber, wave_amplitude=wave_amplitude)

    Qᵘ = EffectiveStressGrowingStokesDrift(wavenumber=wavenumber, wave_amplitude=wave_amplitude,
                                           growth_time_scale=growth_time_scale)

    u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ))

else
    error("There is no such case '$case'!")
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

# Use spinup-save-point command line argument to determine initial condition
iterations = get_iters(filepath)
save_point = args["spinup-save-point"] == 0 ? length(iterations) : args["spinup-save-point"]
set_from_file!(model, filepath, i=save_point)

# # Prepare the simulation

stop_time = 2 * 2π / f

# Adaptive time-stepping
wizard = TimeStepWizard(       cfl = 0.1,
                                Δt = 1e-1,
                        max_change = 1.1,
                            max_Δt = 10.0)

messenger = SimulationProgressMessenger(model, wizard)

simulation = Simulation(model, Δt=wizard, stop_time=stop_time, progress_frequency=100, progress=messenger)

# # Specify output

spinup_iteration = iterations[save_point]
spinup_time = get_time(filepath, spinup_iteration) * f / 2π # Final spinup time / initial sim time, in inertial periods

prefix = @sprintf("%s_Qb%.1e_Nsq%.1e_init%.1f_a%.1f_k%.1e_T%.1f_Nh%d_Nz%d", case,
                  get_parameter(filepath, "boundary_conditions", "Qᵇ"),
                  get_parameter(filepath, "initial_conditions", "N²"),
                  spinup_time, wave_amplitude, wavenumber, 
                  growth_time_scale / hour, model.grid.Nx, model.grid.Nz)

data_directory = joinpath(@__DIR__, "..", "data", prefix) # save data in /data/prefix

# Copy this file into the directory with data
mkpath(data_directory)
cp(@__FILE__, joinpath(data_directory, basename(@__FILE__)), force=true)

"Save a few things that we might want when we analyze the data."
function init(file, model; kwargs...)

    file["sponge_layer/δ"] = δ
    file["sponge_layer/τ"] = τ
    file["initial_conditions/N²"] = N²

    if model.surface_waves !== nothing 
        file["surface_waves/wavenumber"] = wavenumber
        file["surface_waves/wave_amplitude"] = wave_amplitude
        file["surface_waves/growth_time_scale"] = growth_time_scale
    else
        file["surface_waves/wavenumber"] = 0
        file["surface_waves/wave_amplitude"] = 0
        file["surface_waves/growth_time_scale"] = 0
    end

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

# Horizontal averages
averages_writer = JLD2OutputWriter(model, horizontal_averages(model); force=true, init=init,
                                   interval = 10minute,
                                        dir = data_directory,
                                     prefix = prefix * "_averages")

simulation.output_writers[:averages] = averages_writer

# # Run

print_banner(simulation)

run!(simulation)

exit() # Release GPU memory
