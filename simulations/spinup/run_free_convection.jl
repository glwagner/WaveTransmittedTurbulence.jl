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

using Random, Printf, Statistics, ArgParse

"Returns a dictionary of command line arguments."
function parse_command_line_arguments()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--Nh"
            help = "The number of grid points in x, y."
            default = 32

        "--Nz"
            help = "The number of grid points in z."
            default = 32

        "--buoyancy_flux", "-Q"
            help = """The surface buoyancy flux that drives convection in units of m² s⁻³. 
                      A positive buoyancy flux implies cooling."""
            default = 1e-9

        "--buoyancy_gradient", "-Q"
            help = """The buoyancy gradient, or the square of the Brunt-Vaisala frequency N²,
                      at the start of the simulation in units s⁻²."""
            default = 1e-5

        "--device", "-d"
            help = "The CUDA device index on which to run the simulation."
            default = 0
    end

    return parse_args(s)
end

args = parse_command_line_arguments()

@hascuda select_device!(args["device"])

# # Set numerical and physical parameters

# These parameters are set on the command line.
Nh = args["Nh"]                # Number of grid points in x, y
Nz = args["Nz"]                # Number of grid points in z
Qᵇ = args["buoyancy_flux"]     # [m² s⁻³] Buoyancy flux at surface
N² = args["buoyancy_gradient"] # [s⁻²] Initial buoyancy gradient

Lh = 128                   # [m] Grid spacing in x, y (meters)
Lz = 64                    # [m] Grid spacing in z (meters)
θ₀ = 20.0                  # [ᵒC] Surface temperature
 f = 1e-4                  # [s⁻¹] Coriolis parameter

# Create the grid 
grid = RegularCartesianGrid(size=(Nh, Nh, Nz), x=(0, Lh), y=(0, Lh), z=(-Lz, 0))

# Calculate stop time as time when boundary layer depth is h = Lz/2.
# Uses a conservative estimate based on 
#
#   h ∼ √(2 * Qᵇ * stop_time / N²)

stop_time = 4π / f

# # Buoyancy, equation of state, temperature flux, and initial temperature gradient

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(α=2e-4), constant_salinity=35.0)

   Qᶿ = Qᵇ / (buoyancy.gravitational_acceleration * buoyancy.equation_of_state.α)
dθdz₀ = N² / (buoyancy.gravitational_acceleration * buoyancy.equation_of_state.α)
 dθdz = dθdz₀

# # Near-wall LES diffusivity modification + temperature flux specification

# Wall-aware AMD model constant
Δz = Lz/Nz
Cᴬᴹᴰ = SurfaceEnhancedModelConstant(Δz, C₀=1/12, enhancement=2, decay_scale=8Δz)

κₑ_bcs = SurfaceFluxDiffusivityBoundaryConditions(grid, Qᵇ; Cʷ=0.1)

κ₀ = κₑ_bcs.z.top.condition # surface diffusivity
dθdz_surface = - Qᶿ / κ₀    # set temperature gradient = - flux / diffusivity

θ_bcs = TracerBoundaryConditions(grid, top = BoundaryCondition(Gradient, dθdz_surface),
                                       bottom = BoundaryCondition(Gradient, dθdz₀))

# # Sponge layer

δ = 8   # [m] Sponge layer width
τ = 60  # [s] Sponge layer damping time-scale

u_forcing = ParameterizedForcing(Fu, (δ=δ, τ=τ))
v_forcing = ParameterizedForcing(Fv, (δ=δ, τ=τ))
w_forcing = ParameterizedForcing(Fw, (δ=δ, τ=τ))
θ_forcing = ParameterizedForcing(Fθ, (δ=δ, τ=τ, dθdz=dθdz₀, θ₀=θ₀))

# # Model instantiation, initial condition, and model run

prefix = @sprintf("free_convection_Qb%.1e_Nsq%.1e_Nh%d_Nz%d", Qᵇ, N², Nh, Nz)

model = IncompressibleModel(       architecture = has_cuda() ? GPU() : CPU(),
                                           grid = grid,
                                        tracers = :T,
                                       buoyancy = buoyancy,
                                       coriolis = FPlane(f=f),
                                        closure = AnisotropicMinimumDissipation(C=Cᴬᴹᴰ),
                            boundary_conditions = (T=θ_bcs, κₑ=(T=κₑ_bcs,)),
                                        forcing = ModelForcing(u=u_forcing, v=v_forcing, w=w_forcing, T=θ_forcing)
                           )

# Initial condition
ε₀, Δθ, w★ = 1e-6, dθdz₀ * Lz, (Qᵇ * Lz)^(1/3)
Ξ(ε₀, z) = ε₀ * randn() * z / Lz * exp(4z / Lz) # noise
θᵢ(x, y, z) = θ₀ + dθdz₀ * z + Ξ(ε₀ * Δθ, z)
uᵢ(x, y, z) = Ξ(ε₀ * w★, z)

Oceananigans.set!(model, T=θᵢ, u=uᵢ, v=uᵢ, w=uᵢ)

function init(file, model; kwargs...)
    save_global!(file, "sponge_layer", :δ)
    save_global!(file, "sponge_layer", :τ)
    save_global!(file, "initial_conditions", :dθdz)
    save_global!(file, "initial_conditions", :θ₀)
    save_global!(file, "boundary_conditions", :Qᵇ)
    save_global!(file, "boundary_conditions", :Qᶿ)
    return nothing
end

# # Prepare the simulation

# Adaptive time-stepping
wizard = TimeStepWizard(       cfl = 0.2,
                                Δt = 1e-1,
                        max_change = 1.1,
                            max_Δt = 10.0)

messenger = SimulationProgressMessenger(model, wizard)

simulation = Simulation(model, Δt=wizard, stop_time=stop_time, progress_frequency=100, progress=messenger)

# # Output

data_directory = joinpath(@__DIR__, "..", "..", "data", prefix) # save data in /data/prefix

# Three-dimensional field output
fields_to_output = merge(model.velocities, model.tracers, (νₑ=model.diffusivities.νₑ,),
                         prefix_tuple_names(:κₑ, model.diffusivities.κₑ))

field_writer = JLD2OutputWriter(model, FieldOutputs(fields_to_output); force=true, init=init,
                                    interval = (stop_time - wizard.max_Δt) / 10, # output 10 fields
                                max_filesize = 2GiB,
                                         dir = data_directory,
                                      prefix = prefix * "_fields")

# Horizontal averages
averages_writer = JLD2OutputWriter(model, horizontal_averages(model); force=true, init=init,
                                   interval = 1hour, 
                                        dir = data_directory,
                                     prefix = prefix*"_averages")

simulation.output_writers[:fields] = field_writer
simulation.output_writers[:averages] = averages_writer

# # Run

print_banner(simulation)

run!(simulation)

exit() # Release GPU memory
