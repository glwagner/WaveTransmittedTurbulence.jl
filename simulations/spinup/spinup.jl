using Oceananigans, Oceananigans.Diagnostics, Oceananigans.OutputWriters, Random, Printf, Statistics
using Oceananigans: has_cuda

if has_cuda() 
    include("cuda_utils.jl")
    select_device!(1)
end

makeplot = false

include("setup_utils.jl")

Nh, Nz = has_cuda() ?
    (256, 192) :
    (32, 24)

 Lh = 128      # [m] Grid spacing in x, y (meters)
 Lz = 48       # [m] Grid spacing in z (meters)
 N² = 1e-5     # [s⁻²] Initial buoyancy gradient
  f = 1e-4     # [s⁻¹] Coriolis parameter
Qᵘ₀ = 0.0      # [m² s⁻²] Velocity flux / stress at surface
Qᵇ₀ = 1e-9     # [m² s⁻³] Buoyancy flux at surface
  τ = 30       # [s] sponge layer damping time-scale
  δ = 3.0      # [m] sponge layer width
 h₀ = 0.0      # [m] initial mixed layer depth
 tf = 2π/f * 2 # [s] End time for the simulation
 ϵᵘ = 0.0
 ϵᵇ = 1e-6

prefix = @sprintf("spinup_med_convecting_only_Nh%d_Nz%d", Nh, Nz)

@inline ramp(t, δ) = 1 - exp(-t^2 / 2δ^2)
@inline pulse(t, δ) = exp(-t^2 / 2δ^2)

@inline Qᵘ(i, j, grid, time, iter, U, C, p) = p.Qᵘ₀ * ramp(time, p.ramp_up)
@inline Qᵇ(i, j, grid, time, iter, U, C, p) = p.Qᵇ₀ * ramp(time, p.ramp_up)

u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ))
b_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵇ₀), bottom = BoundaryCondition(Gradient, N²))

model = Model(
         architecture = has_cuda() ? GPU() : CPU(),
                 grid = RegularCartesianGrid(size=(Nh, Nh, Nz), length=(Lh, Lh, Lz)),
             buoyancy = BuoyancyTracer(), tracers = :b,
             coriolis = FPlane(f=f),
              closure = AnisotropicMinimumDissipation(),
  boundary_conditions = HorizontallyPeriodicSolutionBCs(b=b_bcs),
              forcing = ModelForcing(u=Fu, v=Fv, w=Fw, b=Fb),
           parameters = (δ=δ, τ=τ, N²=N², h₀=h₀, Qᵘ₀=Qᵘ₀, Qᵇ₀=Qᵇ₀, ramp_up=2π/f)
)

# Initial condition
Ξ(z) = randn() * z / Lz * (1 + z / Lz) # noise
uᵋ(x, y, z) = ϵᵘ * sqrt(abs(Qᵘ₀)) * Ξ(z)
bᵋ(x, y, z) = ϵᵇ * N² * Lz/Nz * Ξ(z)
b₀(x, y, z) = b₀₀(z, N², 0.0) + bᵋ(x, y, z)

Oceananigans.set!(model, u=uᵋ, v=uᵋ, w=uᵋ, b=b₀)
output_writers!(model, prefix)

run_model!(model, tf)

exit() # Release GPU memory
