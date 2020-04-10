function tke(u, v, w)
    U = mean(u, dims=(1, 2)) 
    V = mean(v, dims=(1, 2)) 

    u′ = u .- U
    v′ = v .- V

    u′ = @. @views 0.5 * ( u′[2:end-1, 2:end-1, 2:end-1] + u′[3:end, 2:end-1, 2:end-1] )
    v′ = @. @views 0.5 * ( v′[2:end-1, 2:end-1, 2:end-1] + v′[2:end-1, 3:end, 2:end-1] )
    w′ = @. @views 0.5 * (  w[2:end-1, 2:end-1, 2:end-2] +  w[2:end-1, 2:end-1, 3:end-1] )

    #@show size(u′) size(v′) size(w′)

    e = @. ( u′^2 + v′^2 + w′^2 ) / 2 

    return e
end

function calculate_statistics(grid, u, v, w, b)
    e = tke(u, v, w)

    U = mean(u, dims=(1, 2)) 
    V = mean(v, dims=(1, 2)) 
    B = mean(b, dims=(1, 2)) 

    w² = mean(w.^2, dims=(1, 2)) 
    #w² = 0.5 * (w²[1:end-2] + w²[2:end-1])

    E = dropdims(mean(e, dims=(1, 2)), dims=(1, 2)) 
    b′ = b .- B

    Bz = (B[2:end] - B[1:end-1]) / grid.Δz

    return b′, U, V, B, Bz, E, w², e
end
