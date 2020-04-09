using JLD2

function mixed_layer_depth(B, ΔB, z)
    B .-= @inbounds B[end] # subtract surface value

    i_mld⁺ = searchsortedlast(B, -ΔB)

    Δz = z[2] - z[1]

    if i_mld⁺ < 1
        return 0.0
    else # linearly interpolate to find mixed layer depth
        dBdz = @inbounds (B[i_mld⁺+1] - B[i_mld⁺]) / Δz
        zΔ = @inbounds z[i_mld⁺] - (B[i_mld⁺] + ΔB) / dBdz
        return -zΔ
    end
end

function analyze_averages(filename, uˢ)

    file = jldopen(filename)

    iterations = parse.(Int, keys(file["timeseries/t"]))

    Nz = file["grid/Nz"]
    Δz = file["grid/Δz"]

    nsteps = length(iterations)

    t = zeros(nsteps)

    U = zeros(nsteps, Nz)
    V = zeros(nsteps, Nz)
    B = zeros(nsteps, Nz)

    u² = zeros(nsteps, Nz)
    v² = zeros(nsteps, Nz)
    w² = zeros(nsteps, Nz+1)

    U_uS = zeros(nsteps, Nz)

    ∫U = zeros(nsteps)
    ∫V = zeros(nsteps)

    ∫u² = zeros(nsteps)
    ∫v² = zeros(nsteps)
    ∫w² = zeros(nsteps)

    max_u² = zeros(nsteps)
    max_v² = zeros(nsteps)
    max_w² = zeros(nsteps)

    h = zeros(nsteps)

    for (i, iter) in enumerate(iterations)
        t[i] = file["timeseries/t/$iter"]
        U[i, :] .= dropdims(file["timeseries/U/$iter"], dims=(1, 2))[2:end-1]  
        V[i, :] .= dropdims(file["timeseries/V/$iter"], dims=(1, 2))[2:end-1]
        B[i, :] .= dropdims(file["timeseries/B/$iter"], dims=(1, 2))[2:end-1]
        u²[i, :] .= dropdims(file["timeseries/U²/$iter"], dims=(1, 2))[2:end-1]
        v²[i, :] .= dropdims(file["timeseries/V²/$iter"], dims=(1, 2))[2:end-1]
        w²[i, :] .= dropdims(file["timeseries/W²/$iter"], dims=(1, 2))[2:end]

        U_uS[i, :] .= U[i, :] .* uˢ

        ∫U[i] = sum(U[i, :]) * Δz
        ∫V[i] = sum(V[i, :]) * Δz

        ∫u²[i] = sum(u²[i, :]) * Δz
        ∫v²[i] = sum(v²[i, :]) * Δz
        ∫w²[i] = sum((w²[i, 1:end-1] .+ w²[i, 2:end]) ./ 2) * Δz

        max_u²[i] = maximum(u²[i, :])
        max_v²[i] = maximum(v²[i, :])
        max_w²[i] = maximum(w²[i, :])

        h[i] = mixed_layer_depth(B[i, :], 1e-5, grid.zC)
    end

    return t, U, V, B, U_uS, u², v², w², ∫U, ∫V, ∫u², ∫v², ∫w², max_u², max_v², max_w², h
end

function extra_averages(B, ∫u², ∫v², ∫w²)
    Bz = @. @views (B[:, 2:end] - B[:, 1:end-1]) / Δz
    
    E = @. 1/2 * (∫u² + ∫v² + ∫w²)
    
    return Bz, E
end

function get_averages(filename)

    file = jldopen(filename)

    iterations = parse.(Int, keys(file["timeseries/t"]))

    Nz = file["grid/Nz"]
    Δz = file["grid/Δz"]

    nsteps = length(iterations)

    t = zeros(nsteps)

    U = zeros(nsteps, Nz)
    V = zeros(nsteps, Nz)
    B = zeros(nsteps, Nz)

    u² = zeros(nsteps, Nz)
    v² = zeros(nsteps, Nz)
    w² = zeros(nsteps, Nz+1)
    
    for (i, iter) in enumerate(iterations)
        t[i] = file["timeseries/t/$iter"]
        U[i, :] .= dropdims(file["timeseries/U/$iter"], dims=(1, 2))[2:end-1]  
        V[i, :] .= dropdims(file["timeseries/V/$iter"], dims=(1, 2))[2:end-1]
        B[i, :] .= dropdims(file["timeseries/B/$iter"], dims=(1, 2))[2:end-1]
        u²[i, :] .= dropdims(file["timeseries/U²/$iter"], dims=(1, 2))[2:end-1]
        v²[i, :] .= dropdims(file["timeseries/V²/$iter"], dims=(1, 2))[2:end-1]
        w²[i, :] .= dropdims(file["timeseries/W²/$iter"], dims=(1, 2))[2:end]
    end

    Bz = @. @views (B[:, 2:end] - B[:, 1:end-1]) / Δz

    return t, U, V, B, Bz, u², v², w²
end

function timeseries_from_averages(grid, U, V, B, u², v², w²; mixed_layer_ΔB=1e-5)

    nsteps, Nz = size(U)

    ∫U = zeros(nsteps)
    ∫V = zeros(nsteps)

    ∫u² = zeros(nsteps)
    ∫v² = zeros(nsteps)
    ∫w² = zeros(nsteps)

    h = zeros(nsteps)

    Δz = grid.Δz

    for i = 1:nsteps
        ∫U[i] = sum(U[i, :]) * Δz
        ∫V[i] = sum(V[i, :]) * Δz

        ∫u²[i] = sum(u²[i, :]) * Δz
        ∫v²[i] = sum(v²[i, :]) * Δz
        ∫w²[i] = sum((w²[i, 1:end-1] .+ w²[i, 2:end]) ./ 2) * Δz

        h[i] = mixed_layer_depth(B[i, :], mixed_layer_ΔB, grid.zC)
    end

    return h, ∫U, ∫V, ∫u², ∫v², ∫w²
end

function compute_theoretical_transports(t, U₀, τ)
    Nt = length(t)
    t′ = range(t[1], stop=t[end], length=100Nt)
    Δt′ = t′[2] - t′[1]

    trapz(f, Δx) = 0.5 * sum(f[1:end-1] + f[2:end]) * Δx

    integral = zeros(Complex{Float64}, Nt)

    for i = 2:Nt
        i₁ = searchsortedfirst(t′, t[i-1])
        i₂ = searchsortedfirst(t′, t[i])

        t′′ = t′[i₁:i₂]
        τ′′ = τ.(t′′)

        kernel = @. exp(im*f*t′′) * τ′′

        integral[i] = integral[i-1] + trapz(kernel, Δt′)
    end

    UU = @. U₀ - exp(-im*f*t) * integral
    U = real.(UU)
    V = imag.(UU)

    return U, V
end
