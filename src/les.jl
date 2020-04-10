struct SurfaceEnhancedModelConstant{T} <: Function
             C₀ :: T
             Δz :: T
             z₀ :: T
    enhancement :: T
    decay_scale :: T
end

"""
    SurfaceEnhancedModelConstant(Δz; FT=Float64, C₀=1/12, z₀=-Δz/2,
                                 enhancement=2, decay_scale=8Δz)

Returns a callable object representing a spatially-variable model constant
for an LES eddy diffusivity model with the surface-enhanced form

    ``C(z) = C₀ * (1 + enhancement * exp((z - z₀) / decay_scale)``

"""
function SurfaceEnhancedModelConstant(Δz; FT=Float64, C₀=1/12, z₀=-Δz/2,
                                      enhancement=2, decay_scale=8Δz)

    return SurfaceEnhancedModelConstant{FT}(C₀, Δz, z₀, enhancement, decay_scale)
end

@inline (C::SurfaceEnhancedModelConstant)(x, y, z) =
    C.C₀ * (1 + C.enhancement * exp((z - C.z₀) / C.decay_scale))




