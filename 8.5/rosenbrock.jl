include("newton.jl")
include("bfgs.jl")

r = x -> (1 - x[1])^2 + 100*(x[2] - x[1]^2)^2

grad_r = x -> [
    -2 + 2*x[1] - 400*x[1]*x[2] + 400*x[1]^3;
    200*x[2] - 200*x[1]^2
]

hess_r = x-> [
    2-400*x[2]+1200x[1]^2     -400*x[1];
    -400*x[1]                 200
]

function results(method)
    for m in [0.01, 10, 50]
        x_0 = (rand(2) .+ 1) * m

        ans = newton(r, grad_r, hess_r, x_0)

        println("|x_0| = ", LinearAlgebra.norm(x_0, 2))

        if ans.converged
            println("argument = ", ans.argument, ", iters = ", ans.iters)
        else
            println("not converged")
        end
        println()
    end
end

println("Newton:")
results(newton)

println("\n\nBFGS:")
results(bfgs)