module WaveTransmittedTurbulence

export
    # top level
    select_device!,

    # output.jl
    horizontal_averages,

    # plotting.jl
    make_plot,

    # forcing.jl
    μ,
    Ψᵢ,
    Fu,
    Fv,
    Fw,
    Fθ,
    Fs,
    Fb,

    # setup.jl
    save_global!,
    print_banner,
    prefix_tuple_names,

    # files.jl
    get_iters,
    get_grid,
    get_fields,
    get_wind_stress,
    get_surface_wave_parameters,
    get_parameter,
    set_from_file!,

    # les.jl
    SurfaceEnhancedModelConstant,
    SurfaceFluxDiffusivityBoundaryConditions,
    save_closure_parameters!,

    # progress_messenger.jl
    SimulationProgressMessenger,

    # stokes_drift.jl
    GrowingStokesDrift,
    EffectiveStressGrowingStokesDrift,
    
    # reexport from Oceananigans / CUDAapi
    @hascuda,
    has_cuda

using OffsetArrays, JLD2, Printf

using Oceananigans,
      Oceananigans.AbstractOperations,
      Oceananigans.BoundaryConditions,
      Oceananigans.TurbulenceClosures,
      Oceananigans.Diagnostics,
      Oceananigans.Fields,
      Oceananigans.Operators,
      Oceananigans.Grids,
      Oceananigans.Utils

using Oceananigans.Architectures: device

using Oceananigans.SurfaceWaves: UniformStokesDrift

using Oceananigans.Buoyancy: g_Earth

using Oceananigans: @hascuda, Face, Cell

using GPUifyLoops: @loop, @launch

using CUDAapi: has_cuda

# Don't try to load PyPlot when we don't have a working python.
withplots = false

#=
try
    using PyPlot
    withplots = true
catch; end
=#

macro haspyplot(expr)
    return withplots ? :($(esc(expr))) : :(nothing)
end

@hascuda begin
    using CUDAnative, CUDAdrv

    function select_device!(ndev)
        @show dev = CuDevice(ndev)
        CUDAnative.device!(dev)
        return nothing
    end
end

include("utils.jl")
include("output.jl")
include("forcing.jl")
include("files.jl")
include("les.jl")
include("progress_messenger.jl")

include("stokes_drift.jl")

@haspyplot include("plotting.jl")

end # module
