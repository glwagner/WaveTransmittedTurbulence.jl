using Oceananigans, Oceananigans.Diagnostics, Oceananigans.OutputWriters, Random, Printf, JLD2, Statistics
using Oceananigans.SurfaceWaves: UniformStokesDrift
using Oceananigans: g_Earth, has_cuda

include("setup_utils.jl")

filename = has_cuda() ?
    "data/spinup_convecting_only_Nh256_Nz192_fields_part17.jld2" :
    "data/spinup_Nh32_Nz48_fields.jld2"

file = jldopen(filename)

# Load parameters
Nh = file["grid/Nx"]
Nz = file["grid/Nz"]
Lh = file["grid/Lx"]
Lz = file["grid/Lz"]
Δz = file["grid/Δz"]

 τ = file["sponge_layer/τ"]
 δ = file["sponge_layer/δ"]
h₀ = file["sponge_layer/h₀"]
N² = file["boundary_conditions/N²"]
 f = file["coriolis/f"]

# Load initial condition from file
i = parse(Int, keys(file["timeseries/t"])[end])
 
u₀ = file["timeseries/u/$i"]
v₀ = file["timeseries/v/$i"]
w₀ = file["timeseries/w/$i"]
b₀ = file["timeseries/b/$i"]

close(file)

duration = 2π/f * 2 # [s] Simulation duration in seconds

# Surface waves
const kˢʷ = 2π/100  # [m⁻¹] Surface wave wavenumber
const aˢʷ = 1.5     # [m] Surface wave amplitude
const Tˢʷ = 4hour   # [s] Surface wave growth time scale
const Uˢ = aˢʷ^2 * kˢʷ * sqrt(g_Earth * kˢʷ)
const τˢ = -aˢʷ^2 * sqrt(g_Earth * kˢʷ) / 2Tˢʷ

# Stokes drift vertical and temporal derivative with ramp-up function
@inline ramp(t, δ) = 1 - exp(-t^2 / (2δ^2))
@inline ∂t_ramp(t, δ) = exp(-t^2 / (2δ^2)) * t / δ^2

@inline pulse(t, δ, T) = exp(-(t-T)^2 / (2δ^2))
@inline ∂t_pulse(t, δ, T) = -(t-T) / δ^2 * exp(-(t-T)^2 / (2δ^2))

@inline ∂z_uˢ(z, t) = 2kˢʷ * Uˢ * exp(2kˢʷ * z) * ramp(t, Tˢʷ)
@inline ∂t_uˢ(z, t) =        Uˢ * exp(2kˢʷ * z) * ∂t_ramp(t, Tˢʷ)
@inline ∂z_uˢ_steady(z, t) = 2kˢʷ * Uˢ * exp(2kˢʷ * z)

@inline Qᵘ(i, j, grid, time, iter, U, C, p) = p.τˢ * p.Tˢʷ * ∂t_ramp(time, p.Tˢʷ)

@inline ∂z_uˢ_pulse(z, t) = 2kˢʷ * Uˢ * exp(2kˢʷ * z) * pulse(t, 2Tˢʷ, 8Tˢʷ)
@inline ∂t_uˢ_pulse(z, t) =        Uˢ * exp(2kˢʷ * z) * ∂t_pulse(t, 2Tˢʷ, 8Tˢʷ)
@inline Qᵘ_pulse(i, j, grid, time, iter, U, C, p) = p.τˢ * p.Tˢʷ * ∂t_pulse(time, 2*p.Tˢʷ, 8*Tˢʷ)

function set_initial_condition!(model)

    array_type = typeof(model.velocities.u.data.parent)

    # Set initial condition
    model.velocities.u.data.parent .= array_type(u₀)
    model.velocities.v.data.parent .= array_type(v₀)
    model.velocities.w.data.parent .= array_type(w₀)
    model.tracers.b.data.parent .= array_type(b₀)

    return nothing
end

function wave_forced_model(pre_prefix; pulse=false)

    b_bcs = HorizontallyPeriodicBCs(bottom = BoundaryCondition(Gradient, N²))

    model = Model(
             architecture = has_cuda() ? GPU() : CPU(),
                     grid = RegularCartesianGrid(size=(Nh, Nh, Nz), length=(Lh, Lh, Lz)),
                 buoyancy = BuoyancyTracer(), tracers = :b,
                 coriolis = FPlane(f=f),
                  closure = AnisotropicMinimumDissipation(),
      boundary_conditions = HorizontallyPeriodicSolutionBCs(b=b_bcs),
                  forcing = ModelForcing(u=Fu, v=Fv, w=Fw, b=Fb),
            surface_waves = pulse ? UniformStokesDrift(∂z_uˢ=∂z_uˢ_pulse, ∂t_uˢ=∂t_uˢ_pulse) :
                                    UniformStokesDrift(∂z_uˢ=∂z_uˢ, ∂t_uˢ=∂t_uˢ),
               parameters = (δ=δ, τ=τ, N²=N², h₀=h₀, aˢʷ=aˢʷ, kˢʷ=kˢʷ, Tˢʷ=Tˢʷ)
    )

    if pulse
        prefix = @sprintf("%s_pulse_a%.1f_Nh%d_Nz%d", pre_prefix, aˢʷ, Nh, Nz)
    else
        prefix = @sprintf("%s_ramp_a%.1f_Nh%d_Nz%d", pre_prefix, aˢʷ, Nh, Nz)
    end

    set_initial_condition!(model)

    output_writers!(model, prefix)

    return model
end

function wind_forced_model(pre_prefix; pulse=false)

    if pulse
        u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ_pulse))
    else
        u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ))
    end

    b_bcs = HorizontallyPeriodicBCs(bottom = BoundaryCondition(Gradient, N²))

    model = Model(
             architecture = has_cuda() ? GPU() : CPU(),
                     grid = RegularCartesianGrid(size=(Nh, Nh, Nz), length=(Lh, Lh, Lz)),
                 buoyancy = BuoyancyTracer(), tracers = :b,
                 coriolis = FPlane(f=f),
                  closure = AnisotropicMinimumDissipation(),
      boundary_conditions = HorizontallyPeriodicSolutionBCs(u=u_bcs, b=b_bcs),
                  forcing = ModelForcing(u=Fu, v=Fv, w=Fw, b=Fb),
               parameters = (δ=δ, τ=τ, N²=N², h₀=h₀, aˢʷ=aˢʷ, kˢʷ=kˢʷ, Tˢʷ=Tˢʷ, τˢ=τˢ)
    )

    if pulse
        prefix = @sprintf("%s_pulse_a%.1f_Nh%d_Nz%d", pre_prefix, aˢʷ, Nh, Nz)
    else
        prefix = @sprintf("%s_ramp_a%.1f_Nh%d_Nz%d", pre_prefix, aˢʷ, Nh, Nz)
    end

    set_initial_condition!(model)

    output_writers!(model, prefix)

    return model
end


function wind_with_waves_model(pre_prefix; pulse=false)

    if pulse
        u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ_pulse))
    else
        u_bcs = HorizontallyPeriodicBCs(top = BoundaryCondition(Flux, Qᵘ))
    end

    b_bcs = HorizontallyPeriodicBCs(bottom = BoundaryCondition(Gradient, N²))

    model = Model(
             architecture = has_cuda() ? GPU() : CPU(),
                     grid = RegularCartesianGrid(size=(Nh, Nh, Nz), length=(Lh, Lh, Lz)),
                 buoyancy = BuoyancyTracer(), tracers = :b,
                 coriolis = FPlane(f=f),
                  closure = AnisotropicMinimumDissipation(),
      boundary_conditions = HorizontallyPeriodicSolutionBCs(u=u_bcs, b=b_bcs),
                  forcing = ModelForcing(u=Fu, v=Fv, w=Fw, b=Fb),
            surface_waves = UniformStokesDrift(∂z_uˢ=∂z_uˢ_steady),
               parameters = (δ=δ, τ=τ, N²=N², h₀=h₀, aˢʷ=aˢʷ, kˢʷ=kˢʷ, τˢ=τˢ, Tˢʷ=Tˢʷ)
    )

    if pulse
        prefix = @sprintf("%s_pulse_a%.1f_Nh%d_Nz%d", pre_prefix, aˢʷ, Nh, Nz)
    else
        prefix = @sprintf("%s_ramp_a%.1f_Nh%d_Nz%d", pre_prefix, aˢʷ, Nh, Nz)
    end

    set_initial_condition!(model)

    output_writers!(model, prefix)

    return model
end
