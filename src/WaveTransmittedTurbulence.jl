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
    SurfaceFluxDiffusivityBoundaryConditions,

    # files.jl
    get_iters,
    get_grid,
    get_fields,
    get_wind_stress,
    get_surface_wave_parameters,
    get_parameter,

    # les.jl
    SurfaceEnhancedModelConstant,

    # progress_messenger.jl
    SimulationProgressMessenger,
    
    # reexport
    @hascuda,
    has_cuda

using OffsetArrays, JLD2, Printf

using Oceananigans,
      Oceananigans.AbstractOperations,
      Oceananigans.BoundaryConditions,
      Oceananigans.Diagnostics,
      Oceananigans.Fields,
      Oceananigans.Operators,
      Oceananigans.Grids,
      Oceananigans.Utils

using Oceananigans.Architectures: device

using Oceananigans: @hascuda, Face, Cell

using GPUifyLoops: @loop, @launch

using CUDAapi: has_cuda

function select_device! end

try
    using PyPlot
    withplots = true
catch
    withplots = false
end

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

include("output.jl")
include("forcing.jl")
include("files.jl")
include("setup.jl")
include("les.jl")
include("progress_messenger.jl")

@haspyplot include("plotting.jl")

end # module
