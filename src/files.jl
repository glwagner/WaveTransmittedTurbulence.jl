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

function set_from_file!(model, filename; i=nothing)

    file = jldopen(filename)

    # Load initial condition from file
    iter = i === nothing ? parse(Int, keys(file["timeseries/t"])[end]) :
                           parse(Int, keys(file["timeseries/t"])[i])
    
    u₀ = file["timeseries/u/$iter"]
    v₀ = file["timeseries/v/$iter"]
    w₀ = file["timeseries/w/$iter"]
    b₀ = file["timeseries/b/$iter"]
    
    close(file)

    array_type = typeof(model.velocities.u.data.parent)

    # Set initial condition
    model.velocities.u.data.parent .= array_type(u₀)
    model.velocities.v.data.parent .= array_type(v₀)
    model.velocities.w.data.parent .= array_type(w₀)
    model.tracers.b.data.parent .= array_type(b₀)

end
