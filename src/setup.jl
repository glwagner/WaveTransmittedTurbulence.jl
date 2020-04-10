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
          fields output: $(simulation.output_writers[:fields].filepath)
        averages output: $(simulation.output_writers[:averages].filepath)
    """

    println(banner)

    return nothing
end

function SurfaceFluxDiffusivityBoundaryConditions(grid, Qᵇ; Cʷ=0.1)
    w★ = (Qᵇ * grid.Lz)^(1/3) # surface turbulent velocity scaling
    κ₀ = Cʷ * grid.Δz * w★
    return DiffusivityBoundaryConditions(grid, top = BoundaryCondition(Value, κ₀))
end

function prefix_tuple_names(prefix, tup)
    new_keys = (Symbol(prefix, :_, key) for key in keys(tup))
    return (; zip(new_keys, values(tup))...)
end
