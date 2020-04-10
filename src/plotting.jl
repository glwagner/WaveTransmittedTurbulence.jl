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
