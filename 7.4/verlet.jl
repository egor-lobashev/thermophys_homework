function acceleration(r::Vector)
    return -r / (r[1]^2 + r[2]^2)^1.5
end

function verlet(r_0::Vector, v_0::Vector, tstart::Float64, tend::Float64; nsteps::Integer)
    r = Vector{Vector{Float64}}(undef, nsteps + 1)
    v = Vector{Vector{Float64}}(undef, nsteps + 1)
    r[1] = r_0
    v[1] = v_0
    a = acceleration(r[1])
    a_next = Vector{Float64}(undef, 2)

    trange = range(tstart, tend; length=nsteps+1)
    τ = step(trange)

    for i in 1:nsteps
        r[i+1] = r[i] + τ*v[i] + 0.5*a*τ^2
        a_next = acceleration(r[i+1])
        v[i+1] = v[i] + 0.5*(a + a_next)*τ

        a = a_next
    end

    return trange, [r v]
end