using Oceananigans.AbstractOperations

using Oceananigans: Face, Cell

macro withplots(expr)
    return makeplot ? :($(esc(expr))) : :(nothing)
end

@withplots begin
    using PyPlot

    fig, axs = subplots(ncols=4, figsize=(16, 5))

    """
        makeplot!(axs, model)

    Make a triptych of x-z slices of vertical velocity, temperature, and salinity
    associated with `model` in `axs`.
    """
    function makeplot!(axs, model)
        jhalf = floor(Int, model.grid.Ny/2)

        ## Coordinate arrays for plotting
        xCˣᶻ = repeat(model.grid.xC, 1, model.grid.Nz)
        xFˣᶻ = repeat(model.grid.xF[1:end-1], 1, model.grid.Nz)
        zCˣᶻ = repeat(reshape(model.grid.zC, 1, model.grid.Nz), model.grid.Nx, 1)
        zFˣᶻ = repeat(reshape(model.grid.zF[1:end-1], 1, model.grid.Nz), model.grid.Nx, 1)

        xCˣʸ = repeat(model.grid.xC, 1, model.grid.Ny)
        yCˣʸ = repeat(reshape(model.grid.yC, 1, model.grid.Ny), model.grid.Nx, 1)

        sca(axs[1]); cla()
        title("\$ u \$")
        pcolormesh(xFˣᶻ, zCˣᶻ, Array(interior(model.velocities.u))[:, jhalf, :])
        xlabel("\$ x \$ (m)"); ylabel("\$ z \$ (m)")

        sca(axs[2]); cla()
        title("\$ v \$")
        pcolormesh(xFˣᶻ, zCˣᶻ, Array(interior(model.velocities.v))[:, jhalf, :])
        xlabel("\$ x \$ (m)")

        sca(axs[3]); cla()
        title("\$ w \$")
        pcolormesh(xCˣᶻ, zFˣᶻ, Array(interior(model.velocities.w))[:, jhalf, :])
        xlabel("\$ x \$ (m)")

        sca(axs[4]); cla()
        title("Vertical velocity at \$ z = $(model.grid.zF[Nz-2]) \$ meters")
        pcolormesh(xCˣʸ, yCˣʸ, Array(interior(model.velocities.w))[:, :, Nz-2])
        xlabel("\$ x \$ (m)"); ylabel("\$ y \$ (m)")

        axs[4].yaxis.set_label_position("right")
        axs[4].tick_params(right=true, labelright=true, left=false, labelleft=false)

        [ax.set_aspect(1) for ax in axs]

        tight_layout()

        pause(0.01)

        return nothing
    end
end

function save_if_present!(file, group, model, parameter_name)
    if parameter_name ∈ propertynames(model.parameters)
        file["$group/$parameter_name"] = getproperty(model.parameters, parameter_name)
    end
    return nothing
end

function init_with_bcs(file, model) 

    file["sponge_layer/δ"] = model.parameters.δ
    file["sponge_layer/τ"] = model.parameters.τ
    file["sponge_layer/N²"] = model.parameters.N²
    file["sponge_layer/h₀"] = model.parameters.h₀

    file["boundary_conditions/N²"] = model.boundary_conditions.solution.b.z.bottom.condition

    save_if_present!(file, "boundary_conditions", model, :Qᵇ₀)
    save_if_present!(file, "boundary_conditions", model, :Qᵘ₀)
    save_if_present!(file, "boundary_conditions", model, :Qᵘ₁)
    save_if_present!(file, "boundary_conditions", model, :ramp_up)

    save_if_present!(file, "surface_waves", model, :aˢʷ)
    save_if_present!(file, "surface_waves", model, :kˢʷ)
    save_if_present!(file, "surface_waves", model, :T)

    return nothing
end

function output_writers!(model, prefix)

    f = model.coriolis.f

    if f ≠ 0 
        field_interval = π / (4f)
        averages_interval = π / (100f)
    else
        field_interval = hour
        averages_interval = 5minute
    end

    fields_to_output = merge(model.velocities, model.tracers, 
                             (νₑ=model.diffusivities.νₑ, κₑ=model.diffusivities.κₑ.b))

    field_writer = JLD2OutputWriter(model, FieldOutputs(fields_to_output); interval=field_interval, max_filesize=1GiB,
                                    dir="data", prefix=prefix*"_fields", init=init_with_bcs, force=true)

    # Create scratch space for calculations
    scratch = CellField(model.architecture, model.grid)

    # Extract short field names
    κₑ = model.diffusivities.κₑ.b
    νₑ = model.diffusivities.νₑ
    ph, pn  = model.pressures.pHY′, model.pressures.pNHS
    b = model.tracers.b
    u, v, w = model.velocities

    τ₁₃ = -νₑ * (∂z(u) + ∂x(w))
    τ₂₃ = -νₑ * (∂z(v) + ∂y(w))
    τ₃₃ = -νₑ * ∂z(w)

    ϵʷ = (   νₑ * (∂z(u) + ∂x(w)) * ∂x(w) 
           + νₑ * (∂z(v) + ∂y(w)) * ∂y(w) 
           + 2 * νₑ * ∂z(w)^2 )

    # Define horizontal averages
    U = HorizontalAverage(u)
    V = HorizontalAverage(v)
    B = HorizontalAverage(b)
    NU = HorizontalAverage(νₑ)
    K = HorizontalAverage(κₑ)

    U² = HorizontalAverage(u^2, scratch)
    V² = HorizontalAverage(v^2, scratch)
    W² = HorizontalAverage(w^2, scratch)
    wu = HorizontalAverage(u*w, scratch)
    wv = HorizontalAverage(v*w, scratch)
    wb = HorizontalAverage(b*w, scratch)

    T₁₃ = HorizontalAverage(νₑ * (∂z(u) + ∂x(w)), scratch)
    T₂₃ = HorizontalAverage(νₑ * (∂z(v) + ∂y(w)), scratch)
    q₃ = HorizontalAverage(κₑ * ∂z(b), scratch)

    # Vertical variance terms
      W³ = HorizontalAverage(w^3, scratch)
    wτ₃₃ = HorizontalAverage(-2 * νₑ * w * ∂z(w), scratch)
     w∇p = HorizontalAverage(w * ∂z(ph) + w * ∂z(pn), scratch)
      ϵʷ = HorizontalAverage(ϵʷ, scratch)

      averages = ( 
                   W³ = model ->   W³(model),
                
                    U = model -> U(model),
                    V = model -> V(model),
                    B = model -> B(model),

                   U² = model -> U²(model),
                   V² = model -> V²(model),
                   W² = model -> W²(model),

                   wu = model -> wu(model),
                   wv = model -> wv(model),
                   wb = model -> wb(model),

                   q₃ = model -> q₃(model),
                  τ₁₃ = model -> T₁₃(model),
                  τ₂₃ = model -> T₂₃(model),

                    ν = model -> NU(model),
                    κ = model -> K(model)
               )

    averages_writer = JLD2OutputWriter(model, averages; interval=averages_interval, force=true,
                                       dir="data", prefix=prefix*"_averages", init=init_with_bcs)

    model.output_writers[:fields] = field_writer
    model.output_writers[:averages] = averages_writer

    return averages
end

softstep(x, ξ, δ) = ( tanh( (x-ξ) / δ ) + 1 ) / 2

@inline μ(z, L, δ, τ) = 1/τ * exp(-(z + L)^2 / 2δ^2)
@inline b₀₀(z, N², h) = N² * (z + h)

@inline Fu(i, j, k, grid, time, U, C, p) = @inbounds -μ(grid.zC[k], grid.Lz, p.δ, p.τ) * U.u[i, j, k]
@inline Fv(i, j, k, grid, time, U, C, p) = @inbounds -μ(grid.zC[k], grid.Lz, p.δ, p.τ) * U.v[i, j, k]
@inline Fw(i, j, k, grid, time, U, C, p) = @inbounds -μ(grid.zF[k], grid.Lz, p.δ, p.τ) * U.w[i, j, k]
@inline Fb(i, j, k, grid, time, U, C, p) = 
    @inbounds -μ(grid.zC[k], grid.Lz, p.δ, p.τ) * (C.b[i, j, k] - b₀₀(grid.zC[k], p.N², p.h₀))

function velocity_maxima(model)
    umax = FieldMaximum(abs, model.velocities.u)
    vmax = FieldMaximum(abs, model.velocities.v)
    wmax = FieldMaximum(abs, model.velocities.w)

    return umax, vmax, wmax
end

function print_banner(model)

    banner = """

    Surface wave simulation with

                 Nx, Ny: $(model.grid.Nx)
                     Nz: $(model.grid.Nz)
                 Lx, Ly: $(model.grid.Lx) meters
                     Lz: $(model.grid.Lz) meters
          fields output: $(model.output_writers[:fields].filepath)
        averages output: $(model.output_writers[:averages].filepath)
    """
            
    println(banner)

    return nothing
end

function run_model!(model, end_time; max_Δt=1.0, initial_Δt=0.01, cfl=0.02)

    umax, vmax, wmax = velocity_maxima(model)

    wizard = TimeStepWizard(cfl=cfl, Δt=initial_Δt, max_change=1.1, max_Δt=max_Δt)

    print_banner(model)

    while model.clock.time < end_time

        ## Update the time step associated with `wizard`.
        update_Δt!(wizard, model)

        ## Time step the model forward
        walltime = Base.@elapsed time_step!(model, 100, wizard.Δt)

        ## Print a progress message
        @printf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
                model.clock.iteration, prettytime(model.clock.time), prettytime(wizard.Δt),
                umax(), vmax(), wmax(), prettytime(walltime))

        @withplots makeplot!(axs, model)
    end

    return nothing
end

function fill_solution_halo_regions!(model)
    velocities, tracers = Oceananigans.datatuples(model.velocities, model.tracers)
    Oceananigans.fill_halo_regions!(merge(velocities, tracers), model.boundary_conditions.solution, model.architecture,
                       model.grid, Oceananigans.TimeSteppers.boundary_condition_function_arguments(model)...)
    return nothing
end
