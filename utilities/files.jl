using JLD2, Oceananigans.Grids

function get_iters(filename)
    file = jldopen(filename)
    iters = parse.(Int, keys(file["timeseries/t"]))
    close(file)

    return iters
end

function get_grid(filename)
    file = jldopen(filename)

    Nx = file["grid/Nx"]
    Ny = file["grid/Ny"]
    Nz = file["grid/Nz"]

    Lx = file["grid/Lx"]
    Ly = file["grid/Ly"]
    Lz = file["grid/Lz"]

    close(file)

    grid = RegularCartesianGrid(size=(Nx, Ny, Nz), length=(Lx, Ly, Lz))

    return grid
end

function get_fields(filename, i)
    file = jldopen(filename)

    t = file["timeseries/t/$i"]
    u = file["timeseries/u/$i"]
    v = file["timeseries/v/$i"]
    w = file["timeseries/w/$i"]
    b = file["timeseries/b/$i"]
    νₑ = file["timeseries/νₑ/$i"]
    κₑ = file["timeseries/κₑ/$i"]

    close(file)

    return t, u, v, w, b, νₑ, κₑ
end

function get_wind_stress(filename)
    file = jldopen(filename)
    τ = abs(file["boundary_conditions/Qᵘ₀"])
    close(file)
    return τ
end

function get_surface_wave_parameters(filename)
    file = jldopen(filename)

    aˢʷ = file["surface_waves/aˢʷ"]
    kˢʷ = file["surface_waves/kˢʷ"]

    close(file)

    return aˢʷ, kˢʷ
end

function get_parameter(filename, group, parameter_name)
    parameter = nothing

    jldopen(filename) do file

        if parameter_name ∈ keys(file["$group"])
            parameter = file["$group/$parameter_name"]
        end

    end

    return parameter
end


