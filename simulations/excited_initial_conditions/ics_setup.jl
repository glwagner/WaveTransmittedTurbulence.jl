using Oceananigans, Oceananigans.Diagnostics, Oceananigans.OutputWriters, Random, Printf, Statistics
using Oceananigans.SurfaceWaves: UniformStokesDrift
using Oceananigans: g_Earth, has_cuda

include("setup_utils.jl")

 Lh = 128      # [m] Grid spacing in x, y (meters)
 Lz = 64       # [m] Grid spacing in z (meters)
  f = 1e-4     # [s⁻¹] Coriolis parameter
 Qᵘ = -3.72e-5 # [m² s⁻²] Velocity flux / stress at surface
 Qᵇ = 2.307e-9 # [m³ s⁻²] Buoyancy flux at surface
 N² = 1.936e-5 # [s⁻²] Initial buoyancy gradient
  τ = 10       # [s] sponge layer damping time-scale
  δ = 6.0      # [m] sponge layer width
 h₀ = 0.0      # [m] initial mixed layer depth
 tf = 2π/f * 2 # [s] End time for the simulation
 ϵᵘ = 1e-6
 ϵᵇ = 1e-6

u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ))
b_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵇ),
                                bottom = BoundaryCondition(Gradient, N²))

# Surface waves
const kˢʷ = 0.105   # [m⁻¹] Surface wave wavenumber

const aˢʷ_1x = 0.8              # [m] Surface wave amplitude
const aˢʷ_2x = sqrt(2) * 0.8    # [m] Surface wave amplitude
const aˢʷ_4x = 1.6              # [m] Surface wave amplitude

const Uˢ_1x = aˢʷ_1x^2 * kˢʷ * sqrt(g_Earth * kˢʷ)
const Uˢ_2x = aˢʷ_2x^2 * kˢʷ * sqrt(g_Earth * kˢʷ)
const Uˢ_4x = aˢʷ_4x^2 * kˢʷ * sqrt(g_Earth * kˢʷ)

@inline ∂z_uˢ_1x(z, t) = Uˢ_1x * exp(2kˢʷ * z) * 2kˢʷ
@inline    uˢ_1x(z)    = Uˢ_1x * exp(2kˢʷ * z)

@inline ∂z_uˢ_2x(z, t) = Uˢ_2x * exp(2kˢʷ * z) * 2kˢʷ
@inline    uˢ_2x(z)    = Uˢ_2x * exp(2kˢʷ * z)

@inline ∂z_uˢ_4x(z, t) = Uˢ_4x * exp(2kˢʷ * z) * 2kˢʷ
@inline    uˢ_4x(z)    = Uˢ_4x * exp(2kˢʷ * z)

# Noise
Ξ(z) = randn() * z / Lz * (1 + z / Lz) # noise

uᵋ(x, y, z) = ϵᵘ * sqrt(abs(Qᵘ)) * Ξ(z)
bᵋ(x, y, z) = ϵᵇ * N² * Lz/Nz * Ξ(z)
b₀(x, y, z) = b₀₀(z, N², h₀) + bᵋ(x, y, z)

function initial_conditions_model(experiment, aˢʷ, Nh, Nz)

    if aˢʷ == 0.8

        ∂z_uˢ = ∂z_uˢ_1x
           uˢ = uˢ_1x
           Uˢ = Uˢ_1x

        prefix = experiment == "resting" ? 
            @sprintf("ics_resting_1x_Nh%d_Nz%d", Nh, Nz) : 
            @sprintf("ics_excited_1x_Nh%d_Nz%d", Nh, Nz)

    elseif aˢʷ == sqrt(2) * 0.8

        ∂z_uˢ = ∂z_uˢ_2x
           uˢ = uˢ_2x
           Uˢ = Uˢ_2x

        prefix = experiment == "resting" ? 
            @sprintf("ics_resting_2x_Nh%d_Nz%d", Nh, Nz) :
            @sprintf("ics_excited_2x_Nh%d_Nz%d", Nh, Nz)

    elseif aˢʷ == 1.6

        ∂z_uˢ = ∂z_uˢ_4x
           uˢ = uˢ_4x
           Uˢ = Uˢ_4x

        prefix = experiment == "resting" ? 
            @sprintf("ics_resting_4x_Nh%d_Nz%d", Nh, Nz) :
            @sprintf("ics_excited_4x_Nh%d_Nz%d", Nh, Nz)

    end 

    if experiment == "control"
        prefix = @sprintf("ics_control_Nh%d_Nz%d", Nh, Nz)
        Uˢ = 0.0
    end

    @show experiment Uˢ

    model_setup = (architecture = has_cuda() ? GPU() : CPU(),
                           grid = RegularCartesianGrid(size=(Nh, Nh, Nz), length=(Lh, Lh, Lz)),
                       buoyancy = BuoyancyTracer(), tracers = :b,
                       coriolis = FPlane(f=f),
                        closure = AnisotropicMinimumDissipation(),
                        forcing = ModelForcing(u=Fu, v=Fv, w=Fw, b=Fb),
            boundary_conditions = HorizontallyPeriodicSolutionBCs(u=u_bcs, b=b_bcs))
    
    if experiment == "control"
        model = Model(; model_setup...,
                         parameters = (δ=δ, τ=τ, N²=N², h₀=h₀, Qᵘ₀=Qᵘ, aˢʷ=0, kˢʷ=0))
    else
        model = Model(; model_setup...,
                      surface_waves = UniformStokesDrift(∂z_uˢ=∂z_uˢ),
                         parameters = (δ=δ, τ=τ, N²=N², h₀=h₀, Qᵘ₀=Qᵘ, aˢʷ=aˢʷ, kˢʷ=kˢʷ))
    end

    u₀(x, y, z) = experiment == "excited" ? uˢ(z) + uᵋ(x, y, z) : uᵋ(x, y, z)

    Oceananigans.set!(model, u=u₀, v=uᵋ, w=uᵋ, b=b₀)

    output_writers!(model, prefix)

    return model
end

