using PyPlot, PyCall, Statistics

mplot3d = pyimport("mpl_toolkits.mplot3d")
Axes3D = mplot3d.Axes3D
matplotlib = pyimport("matplotlib")
ticker = pyimport("matplotlib.ticker")
ScalarFormatter = ticker.ScalarFormatter

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

function tke(u, v, w)
    U = mean(u, dims=(1, 2))
    V = mean(v, dims=(1, 2))

    u′ = u .- U
    v′ = v .- V

    u′ = @. @views 0.5 * ( u′[2:end-1, 2:end-1, 2:end-1] + u′[3:end, 2:end-1, 2:end-1] )
    v′ = @. @views 0.5 * ( v′[2:end-1, 2:end-1, 2:end-1] + v′[2:end-1, 3:end, 2:end-1] )
    w′ = @. @views 0.5 * (  w[2:end-1, 2:end-1, 2:end-1] +  w[2:end-1, 2:end-1, 3:end] )

    e = @. ( u′^2 + v′^2 + w′^2 ) / 2

    return e
end

function calculate_statistics(u, v, w, b)
    e = tke(u, v, w)

    U = mean(u, dims=(1, 2))
    V = mean(v, dims=(1, 2))
    B = mean(b, dims=(1, 2))

    w² = mean(w.^2, dims=(1, 2))
    w² = 0.5 * (w²[1:end-2] + w²[2:end-1])

    E = dropdims(mean(e, dims=(1, 2)), dims=(1, 2))
    b′ = b .- B

    Bz = (B[2:end] - B[1:end-1]) / grid.Δz

    return b′, U, V, B, Bz, E, w², e
end

get_position(ax) = [b for b in ax.get_position().bounds]
