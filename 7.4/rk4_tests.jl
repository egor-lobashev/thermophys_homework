include("rk4.jl")

function func(t, u)
    return [
        u[1]+u[2];
        u[1]+u[2]+t
    ]
end

problem = CauchyODEProblem(; f=func, tstart=0, tend=1, u₀=[0; 0])

steps = 10
trange, u = rk4(problem, nsteps=steps)
for i in 1:(steps+1)
    println(trange[i], ", ", u[i][1], ", ", u[i][2])
end

println()
###############################

function another_func(t, u)
    return [
        u[2]
        -u[2]+2t
    ]
end

another_problem = CauchyODEProblem(; f=another_func, tstart=0, tend=1, u₀=[0; 0])

steps = 10
trange, u = rk4(another_problem, nsteps=steps)
for i in 1:(steps+1)
    println(trange[i], ", ", u[i][1])
end