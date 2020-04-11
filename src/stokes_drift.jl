uˢ(a, k, g=g_Earth) = a^2 * k * √(g * k)

struct GrowingStokesShear{T} <: Function
    wavenumber :: T
    wave_amplitude :: T
    max_surface_drift :: T
    growth_time_scale:: T
end

struct GrowingStokesTendency{T} <: Function
    wavenumber :: T
    wave_amplitude :: T
    max_surface_drift :: T
    growth_time_scale :: T
end

# Stokes drift vertical and temporal derivative with ramp-up function
@inline ramp(t, T) = 1 - exp(-t^2 / (2T^2))
@inline ∂t_ramp(t, T) = exp(-t^2 / (2T^2)) * t / T^2

@inline ∂z_uˢ_growing(z, t, k, U, T) = 2k * U * exp(2k * z) * ramp(t, T)
@inline ∂t_uˢ_growing(z, t, k, U, T) =      U * exp(2k * z) * ∂t_ramp(t, T)


@inline (uˢ::GrowingStokesShear)(z, t) = ∂z_uˢ_growing(z, t, uˢ.wavenumber, uˢ.max_surface_drift,
                                                       uˢ.growth_time_scale)

@inline (uˢ::GrowingStokesTendency)(z, t) = ∂t_uˢ_growing(z, t, uˢ.wavenumber, uˢ.max_surface_drift,
                                                          uˢ.growth_time_scale)

growing_wind_stress_func(x, y, t, p) = p.τʷ * p.Tʷ * ∂t_ramp(t, p.Tʷ)

function GrowingStokesDrift(; wavenumber, wave_amplitude, growth_time_scale)

    # Deep water surface wave Stokes drift amplitude...
    max_surface_drift = uˢ(wave_amplitude, wavenumber)

    shear = GrowingStokesShear(wavenumber, wave_amplitude, max_surface_drift, growth_time_scale)
    tendency = GrowingStokesTendency(wavenumber, wave_amplitude, max_surface_drift, growth_time_scale)

    return UniformStokesDrift(∂z_uˢ=shear, ∂t_uˢ=tendency)
end

#####
##### Effective stress
#####

function EffectiveStressGrowingStokesDrift(; wavenumber, wave_amplitude, growth_time_scale)

    effective_stress = - wave_amplitude^2 * sqrt(g_Earth * wavenumber) / (2 * growth_time_scale)

    growing_wind_stress = BoundaryFunction{:z, Face, Cell}(growing_wind_stress_func,
                                                           (τʷ = effective_stress, Tʷ = growth_time_scale))

    return growing_wind_stress
end

#####
##### Steady Stokes shear
#####

struct SteadyStokesShear{T} <: Function
    wavenumber :: T
    wave_amplitude :: T
    surface_drift :: T
end

@inline ∂z_uˢ_steady(z, t, k, U) = 2k * U * exp(2k * z)

@inline (uˢ::SteadyStokesShear)(z, t) = ∂z_uˢ_steady(z, t, uˢ.wavenumber, uˢ.surface_drift)

function SteadyStokesDrift(; wavenumber, wave_amplitude)

    # Deep water surface wave Stokes drift amplitude...
    surface_drift = wave_amplitude^2 * wavenumber * sqrt(g_Earth * wavenumber)

    shear = SteadyStokesShear(wavenumber, wave_amplitude, surface_drift)

    return UniformStokesDrift(∂z_uˢ=shear)
end
