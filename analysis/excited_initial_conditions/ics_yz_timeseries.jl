include("visualization.jl")
include("files.jl")

using Printf

suffix = "4x"

r1_filename = "data/ics_resting_4x_Nh256_Nz256_fields_part2.jld2"
e1_filename = "data/ics_excited_4x_Nh256_Nz256_fields_part2.jld2"

f = get_parameter(r1_filename, "coriolis", "f")
τ = get_wind_stress(r1_filename)

grid = get_grid(r1_filename)

xCˣᶻ = repeat(grid.xC, 1, grid.Nz)
xFˣᶻ = repeat(grid.xF[1:end-1], 1, grid.Nz)
zCˣᶻ = repeat(reshape(grid.zC, 1, grid.Nz), grid.Nx, 1)
zFˣᶻ = repeat(reshape(grid.zF[1:end-1], 1, grid.Nz), grid.Nx, 1)

fs = 10
hpad = 0.02
vpad = 0.06
zlim = -32
hlim = 128
box_props = Dict(:alpha=>0.8, :facecolor=>"w", :edgecolor=>"w")

wlim = 5
wlevels = vcat([-2wlim], -wlim:2wlim/11:wlim, [2wlim])

i = round(Int, grid.Nx/2)
j = round(Int, grid.Ny/2)

function plot_w_contours!(ax1, ax2, filename, name; n=nothing)
    iters = get_iters(filename)
    iter = n === nothing ? iters[end] : iters[n]
    t, u, v, w, b, ν, κ = get_fields(filename, iter)

    @show f * t / 2π

    sca(ax1)
    image = contourf(xCˣᶻ, zFˣᶻ, w[2:end-1, j, 2:end-1] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

    text(hpad*grid.Lx, zlim+vpad*grid.Lz, name, fontsize=fs, bbox=box_props)
    xlim(0, hlim)
    ylim(zlim, -grid.Δz)

    sca(ax2)
    contourf(xCˣᶻ, zFˣᶻ, w[i, 2:end-1, 2:end-1] ./ sqrt(τ), cmap="RdBu_r", levels=wlevels, vmin=-wlim, vmax=wlim)

    text(hpad*grid.Lx, zlim+vpad*grid.Lz, name, fontsize=fs, bbox=box_props)
    xlim(0, hlim)
    ylim(zlim, -grid.Δz)

    return image
end

close("all")

fig, axs = subplots(ncols=2, nrows=2, figsize=(9, 6))

im_r = plot_w_contours!(axs[1, 1], axs[2, 1], r1_filename, "Equilibrated", n=1)
im_e = plot_w_contours!(axs[1, 2], axs[2, 2], e1_filename, "Excited", n=1)

ni, nj = size(axs)

for i = 1:ni
    sca(axs[i, 1])
    ylabel("\$ z \$ (m)")

    axs[i, 2].tick_params(left=false, labelleft=false)
end

for j = 1:nj
    sca(axs[1, j])
    xlabel("\$ x \$ (m)")

    sca(axs[2, j])
    xlabel("\$ y \$ (m)")

end

for ax in axs
    ax.set_aspect(1)
end

get_bounds(ax) = [b for b in ax.get_position().bounds]

pos = Array{Array{Float64, 1}, 2}(undef, ni, nj)

#=
for i = 1:3, j=1:4
    pos[i, j] = get_bounds(axs[i, j])
end

Δx = pos[1, 2][1] - pos[1, 1][1]
Δy = pos[1, 1][2] - pos[2, 1][2]

yshift!(pos, shift) = pos[2] += shift * Δy
xshift!(pos, shift) = pos[1] += shift * Δx

yshift = 0.1
xshift = 0.1

for j=1:4
    yshift!(pos[1, j], -yshift)
    yshift!(pos[3, j], yshift)
end

for i = 1:3
    xshift!(pos[i, 2], -xshift)
    xshift!(pos[i, 3], -xshift)
    xshift!(pos[i, 4], -2xshift)
end

for i = 1:3, j=1:4
    axs[i, j].set_position(pos[i, j])
end

colorbar(image, ax=axs, ticks=-wlim:2wlim/5:wlim, pad=0.01)

savefig(@sprintf("w_xz_yz_timeseries_%s.png", suffix), dpi=480)
=#
