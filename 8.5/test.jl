include("newton.jl")
include("bfgs.jl")

function results(method)
    for n in [2, 10, 100]
        B = rand(n, n)
        A = B*B'

        println("---------------------------\ndimensions = ", n)

        for m in [0.01, 10, 1000, 1000000]
            x_0 = (rand(n) .+ 1) * m

            ans = newton(x->dot(x, A*x), x->2*A*x, x->2*A, x_0)

            println("|x_0| = ", LinearAlgebra.norm(x_0, 2))

            if ans.converged
                println("|argument| = ", LinearAlgebra.norm(ans.argument), ", iters = ", ans.iters)
            else
                println("not converged")
            end
            println()
        end
    end
end

println("Newton:")
results(newton)

println("\n\nBFGS:")
results(bfgs)