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
        "--Nh"
            help = "The number of grid points in x, y."
            default = 32
            arg_type = Int

        "--Nz"
            help = "The number of grid points in z."
            default = 32
            arg_type = Int

        "--initial_condition"
            help = """The initial conditions to use. The options are
                   1. "excited"
                   2. "resting"
                   """
            default = "excited"
            arg_type = String

        "--wave_multiplier"
            help = "Sets the wave height as wave_multiplier * 0.8 meters."
            default = 1.0
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
 
# # Specify numerical and physical parameters

Nh = args["Nh"]      # Number of grid points in x, y
Nz = args["Nz"]      # Number of grid points in z

# Coriolis parameter, kinematic stress, buoyancy flux, and wave parameters from 
#
#   > McWilliams et al, "Langmuir Turbulence in the Ocean," JFM (1997)
#
# Note that our domain is not the same size as McWilliams et al. (1997).

        Lh = 128       # [m] Grid spacing in x, y (meters)
        Lz = 64        # [m] Grid spacing in z (meters)
         f = 1e-4      # [s⁻¹] Coriolis parameter
        Qᵘ = -3.72e-5  # [m² s⁻²] Velocity flux / stress at surface
        Qᵇ = 2.307e-9  # [m³ s⁻²] Buoyancy flux at surface
        N² = 1.936e-5  # [s⁻²] Initial buoyancy gradient
wavenumber = 0.105     # [m⁻¹] Wavenumber of the steady monchromatic wave field overhead
 stop_time = 4π / f    # [s] End time for the simulation

# A 'wave multiplier' of 1 corresponds to original parameters from McWilliams et al (1997).
# A multipler of 0 means no waves, and > 1 produces strog waves
wave_amplitude = 0.8 * args["wave_multiplier"]

# # Choose initial condition

ic = args["initial_condition"]

if ic == "resting"
    uᵢ(x, y, z) = 0.0
elseif ic == "excited"
    uᵢ(x, y, z) = uˢ(wave_amplitude, wavenumber) * exp(2 * wavenumber * z)
end

# # Choose Stokes drift

if wave_amplitude == 0
    stokes_drift = nothing # minor optimziation for control run.
else
    stokes_drift = SteadyStokesDrift(wavenumber=wavenumber, wave_amplitude=wave_amplitude)
end

# # Set up sponge layer

τ = 60       # [s] sponge layer damping time-scale
δ = 8        # [m] sponge layer width

u_forcing = ParameterizedForcing(Fu, (δ=δ, τ=τ))
v_forcing = ParameterizedForcing(Fv, (δ=δ, τ=τ))
w_forcing = ParameterizedForcing(Fw, (δ=δ, τ=τ))
b_forcing = ParameterizedForcing(Fb, (δ=δ, τ=τ, dbdz=N²))

# # Set up boundary conditions

grid = RegularCartesianGrid(size=(Nh, Nh, Nz), length=(Lh, Lh, Lz))

u_bcs = UVelocityBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵘ))
b_bcs = TracerBoundaryConditions(grid, top = BoundaryCondition(Flux, Qᵇ))

# # Model instantiation and initial condition

model = IncompressibleModel(       architecture = has_cuda() ? GPU() : CPU(),
                                           grid = grid,
                                        tracers = :b,
                                       buoyancy = BuoyancyTracer(),
                                       coriolis = FPlane(f=f),
                                  surface_waves = stokes_drift,
                                        closure = AnisotropicMinimumDissipation(),
                            boundary_conditions = (b=b_bcs, u=u_bcs),
                                        forcing = ModelForcing(u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)
                           )


# Noise
Ξ(z) = randn() * z / Lz * (1 + z / Lz) # noise
ϵᵘ, ϵᵇ = 1e-6, 1e-6

uᵋ(x, y, z) = ϵᵘ * sqrt(abs(Qᵘ)) * Ξ(z)
bᵋ(x, y, z) = ϵᵇ * N² * Lz / Nz * Ξ(z)
bᵢ(x, y, z) = N² * z + bᵋ(x, y, z)

set!(model, v=uᵋ, w=uᵋ, b=bᵢ, u = (x, y, z) -> uᵢ(x, y, z) + uᵋ(x, y, z))

# # Prepare the simulation

stop_time = 4π / f

# Adaptive time-stepping
wizard = TimeStepWizard(       cfl = 0.1,
                                Δt = 1e-1,
                        max_change = 1.1,
                            max_Δt = 10.0)

messenger = SimulationProgressMessenger(model, wizard)

simulation = Simulation(model, Δt=wizard, stop_time=stop_time, progress_frequency=100, progress=messenger)

# # Specify output

prefix = @sprintf("initial_condition_study_%s_%sx_Nh%d_Nz%d", ic, string(args["wave_multiplier"]), Nh, Nz)

data_directory = joinpath(@__DIR__, "..", "data", prefix) # save data in /data/prefix

# Copy this file into the directory with data
mkpath(data_directory)
cp(@__FILE__, joinpath(data_directory, basename(@__FILE__)), force=true)

"Save a few things that we might want when we analyze the data."
function init(file, model; kwargs...)
    file["sponge_layer/δ"] = δ
    file["sponge_layer/τ"] = τ
    file["initial_conditions/N²"] = N²
    file["boundary_conditions/Qᵇ"] = Qᵇ
    file["boundary_conditions/Qᵘ"] = Qᵘ
    file["surface_waves/wavenumber"] = wavenumber
    file["surface_waves/wave_amplitude"] = wave_amplitude
    return nothing
end

# Three-dimensional field output
fields_to_output = merge(model.velocities, model.tracers, (νₑ=model.diffusivities.νₑ,),
                         prefix_tuple_names(:κₑ, model.diffusivities.κₑ))

field_writer = JLD2OutputWriter(model, FieldOutputs(fields_to_output); force=true, init=init,
                                    interval = π / 2f,
                                max_filesize = 2GiB,
                                         dir = data_directory,
                                      prefix = prefix * "_fields")

simulation.output_writers[:fields] = field_writer

# Horizontal averages
#=
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
