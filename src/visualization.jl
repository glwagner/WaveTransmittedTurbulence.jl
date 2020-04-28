function bottom_left_text!(ax, txt; kwargs...)
    sca(ax)
    text(0.015, 0.02, txt; transform=ax.transAxes, horizontalalignment="left", verticalalignment="bottom", kwargs...)
    return nothing
end

function sequential_limit_levels(field; saturate=0.6, nlevs=15, limit=nothing)
    max_field = maximum(field)
    limit = limit === nothing ? saturate * max_field : limit

    levels = max_field > limit ?
                levels = vcat(range(0.0, stop=limit, length=nlevs), [max_field]) :
                levels = vcat(range(0.0, stop=limit, length=nlevs))

    return limit, levels
end

function divergent_limit_levels(field; saturate=0.6, nlevs=14, limit=nothing)
    max_field = maximum(abs, field)
    limit = limit === nothing ? saturate * max_field : limit

    levels = max_field > limit ?
                levels = vcat([-max_field], range(-limit, stop=limit, length=nlevs), [max_field]) :
                levels = range(-limit, stop=limit, length=nlevs)

    return limit, levels
end

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

function usecmbright()
    rc("text.latex", preamble="\\usepackage{cmbright}")
    rc("font", family="sans-serif")
    return nothing
end

defaultcolors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

function removespine(ax, side)
    ax.spines[side].set_visible(false)
    return nothing
end

removespines(ax::PyObject, sides...) = [removespine(ax, side) for side in sides]
removespines(sides...) = [removespine(gca(), side) for side in sides]

axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
make_axes_locatable = axes_grid1.make_axes_locatable

function meshgrid(x, y)
    X = repeat(x, 1, length(y))
    Y = repeat(reshape(y, 1, length(y)), length(x), 1)
    return X, Y
end

get_position(ax) = [b for b in ax.get_position().bounds]

function plot_profiles!(axs, grid, B, Bz, S, w²; label="", kwargs...)
    if B != nothing
        sca(axs[1])
        plot(B, grid.zC; label=label, kwargs...)
    end

    if Bz != nothing
        sca(axs[2])
        plot(Bz[2:end-1], grid.zF[2:end-1]; kwargs...)
    end

    if S != nothing
        sca(axs[3])
        plot(S,  grid.zC; kwargs...)
    end

    if w² != nothing
        sca(axs[4])
        plot(w², grid.zF; kwargs...)
    end

    return nothing
end

function shift_up!(ax, shift)
    pos = get_position(ax)
    pos[2] += shift
    ax.set_position(pos)
    return nothing
end

function shift_down!(ax, shift)
    pos = get_position(ax)
    pos[2] -= shift
    ax.set_position(pos)
    return nothing
end

function shift_right!(ax, shift)
    pos = get_position(ax)
    pos[1] += shift
    ax.set_position(pos)
    return nothing
end

function shift_left!(ax, shift)
    pos = get_position(ax)
    pos[1] -= shift
    ax.set_position(pos)
    return nothing
end

function stretch_x!(ax, stretch)
    pos = get_position(ax)
    pos[3] += stretch
    ax.set_position(pos)
    return nothing
end

function stretch_y!(ax, stretch)
    pos = get_position(ax)
    pos[4] += stretch
    ax.set_position(pos)
    return nothing
end
