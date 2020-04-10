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
