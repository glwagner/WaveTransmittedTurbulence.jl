using Oceananigans, Printf

include("analysis.jl")
include("visualization.jl")
include("files.jl")

suffix = "1x"
Nh, Nz = 256, 256
control_filename = @sprintf("data/ics_control_Nh%d_Nz%d_averages.jld2", Nh, Nz)
resting_filename = @sprintf("data/ics_resting_%s_Nh%d_Nz%d_averages.jld2", suffix, Nh, Nz) 
excited_filename = @sprintf("data/ics_excited_%s_Nh%d_Nz%d_averages.jld2", suffix, Nh, Nz)

grid = get_grid(resting_filename)

t_c, U_c, V_c, B_c, Bz_c, u²_c, v²_c, w²_c = get_averages(control_filename)
t_r, U_r, V_r, B_r, Bz_r, u²_r, v²_r, w²_r = get_averages(resting_filename)
t_e, U_e, V_e, B_e, Bz_e, u²_e, v²_e, w²_e = get_averages(excited_filename)

#####
##### Plot
#####

close("all")
fig, axs = subplots(ncols=4, figsize=(10, 4), sharey=true)

α₁ = 0.6
α₂ = 1.0

lw₂ = 2.5
lw₁ = 1.5

depth = 36

k = searchsortedfirst(grid.zF, -depth)

@show length(t_c)
@show length(t_r)
@show length(t_e)

iters = [51, 401]

c_labels, r_labels, e_labels = Tuple([] for i = 1:3)

f = get_parameter(resting_filename, "coriolis", "f")

for (j, i) in enumerate(iters)
    @show t_c[i] / hour t_c[i] * 1e-4 / 2π

    if j == 1
        if isapprox(f * t_c[i] / 2π, 0.25, atol=1e-2)
            push!(c_labels, @sprintf("No waves, \$ t = \\frac{1}{4} \\times 2\\pi/f\$"))
        else
            push!(c_labels, @sprintf("No waves, \$ t = %.2f \\times 2\\pi/f\$", f * t_c[i] / 2π))
        end
    else
        if isapprox(f * t_c[i] / 2π, 2.0, atol=1e-2)
            push!(c_labels, @sprintf("\$ t = %d \\times 2\\pi/f\$", f * t_c[i] / 2π))
        else
            push!(c_labels, @sprintf("\$ t = %.2f \\times 2\\pi/f\$", f * t_c[i] / 2π))
        end
    end

    push!(r_labels, "Equilibrated")
    push!(e_labels, "Excited")
end

for (j, i) in enumerate(iters)

    r_label = r_labels[j]
    e_label = e_labels[j]

    sca(axs[1])
    plot(B_c[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₁, "-")  
    plot(B_r[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, "--")
    plot(B_e[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, ":") 

    sca(axs[2])
    plot(Bz_c[i, :], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₂, alpha=α₁, "-")
    plot(Bz_r[i, :], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₂, "--")
    plot(Bz_e[i, :], grid.zF[2:end-1], color=defaultcolors[j], linewidth=lw₁, alpha=α₂, ":")

    sca(axs[3])
    plot(U_c[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₁, label=c_labels[j], "-")
    plot(U_r[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, label=r_label, "--")
    plot(U_e[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, label=e_label, ":")

    sca(axs[4])
    plot(V_c[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₂, alpha=α₁, "-")
    plot(V_r[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, "--")
    plot(V_e[i, :], grid.zC, color=defaultcolors[j], linewidth=lw₁, alpha=α₂, ":")

end

sca(axs[1])
xlabel(L"\bar b \, (\mathrm{m \, s^{-2}})")
ylabel("\$ z \$ (m)")

sca(axs[2])
xlabel(L"\partial_z \bar b \, (\mathrm{s^{-2}})")

sca(axs[3])
xlabel(L"\bar u \, (\mathrm{m \, s^{-1}})")
legend(loc=3, prop=Dict("size" => 12), bbox_to_anchor=(0.3, 0., 1.0, 0.5))

sca(axs[4])
xlabel(L"\bar v \, (\mathrm{m \, s^{-1}})")
ylabel("\$ z \$ (m)")

removespines(axs[1], "top", "right")
removespines(axs[end], "top", "left")

for j = 2:3
    removespines(axs[j], "top", "right", "left")
    axs[j].tick_params(left=false, labelleft=false)
end

axs[end].yaxis.set_label_position("right")
axs[end].tick_params(left=false, labelleft=false, right=true, labelright=true)
axs[end].set_zorder(-1)

positions = [get_position(ax) for ax in axs]

xwiggle = 0.01
positions[1][1] -= 2xwiggle
positions[2][1] -= xwiggle
positions[3][1] += xwiggle
positions[4][1] += 2xwiggle

xshift = 0.015
yshift = 0.05
for (i, pos) in enumerate(positions)
    pos[1] -= xshift
    pos[2] += yshift
    axs[i].set_position(pos)
end

sca(axs[1])
xlim(B_c[1, k], B_c[1, end])
ylim(-depth, 0.1)

#axs[1].xaxis.set_ticks([-5e-4, -2e-4]) #set_major_formatter(formatter)

formatter = ScalarFormatter() #useMathText=true)
formatter.set_powerlimits((-2, 2))
axs[1].xaxis.set_major_formatter(formatter)

formatter = ScalarFormatter() #useMathText=true)
formatter.set_powerlimits((-2, 2))
axs[2].xaxis.set_major_formatter(formatter)


savefig("ics_averages_$suffix.png", dpi=480)
