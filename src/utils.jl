function save_global!(file, group, name)
    val = try
        eval(name)
    catch err
        elaboration = err isa UndefVarError ? "because $(err.var) is not defined in the global scope." : ""
        @warn "$name will not be saved $elaboration"
        nothing
    end

    if val !== nothing
        file["$group/$name"] = val
    end

    return nothing
end

function print_banner(simulation)

    model = simulation.model

    banner = """

    Simulation with

                 Nx, Ny: $(model.grid.Nx)
                     Nz: $(model.grid.Nz)
                 Lx, Ly: $(model.grid.Lx) meters
                     Lz: $(model.grid.Lz) meters
                 output: $(simulation.output_writers[:fields].filepath)
    """
        #averages output: $(simulation.output_writers[:averages].filepath)

    println(banner)

    return nothing
end

function prefix_tuple_names(prefix, tup)
    new_keys = (Symbol(prefix, :_, key) for key in keys(tup))
    return (; zip(new_keys, values(tup))...)
end
