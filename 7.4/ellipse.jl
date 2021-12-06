include("rk4.jl")
include("verlet.jl")

function func(t, u)
    r3 = (u[1]^2 + u[2]^2)^1.5
    return [
        u[3]
        u[4]
        -u[1]/r3
        -u[2]/r3
    ]
end

E = u -> (u[3]^2 + u[4]^2)/2 - (u[1]^2 + u[2]^2)^-0.5

r_0 = -1.0
v_0 = 1.2
u_0 = [r_0; 0; 0; v_0]

E_0 = E(u_0)
# println(E_0)
a = -1/(2*E_0)
T = 2*pi*a^1.5

t_end = T * 3
steps = 20 * 3

problem = CauchyODEProblem(; f=func, tstart=0, tend=t_end, u₀=u_0)
trange_rk, u_rk = rk4(problem, nsteps=steps)
trange_v, u_v = verlet([r_0; 0], [0; v_0], 0.0, t_end, nsteps=steps)

function print_results(trange, u, verlet=false)
    for i in 1:(steps+1)
        ui = Vector{Float64}(undef, 4)
        if verlet
            ui = [u[i, 1]; u[i, 2]]
        else
            ui = u[i]
        end
        
        println(u[i][1], ", ", u[i][2])
        # println(trange[i], ", ", E(ui))
    end
end

print_results(trange_rk, u_rk)
println()
print_results(trange_v, u_v, true)