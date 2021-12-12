include("cholesky.jl")
include("backtracking.jl")

struct BFGSResult{T<:Real}
    converged::Bool      # метод сошёлся
    argument::Vector{T}  # найденный минимум
    iters::Int           # число прошедших итераций
end

function bfgs(f, ∇f, x0, invH0=diagm(0=>ones(length(x)));
    maxiter=200,
    gtol=1e-5,
)
    x = float.(x0)
    
    # проверяем сходимость по норме-2 градиента
    g = ∇f(x)

    if LinearAlgebra.norm(g, 2) <= gtol
        return BFGSResult(true, x, 0)
    end

    # обратный квази-гессиан
    invB = Matrix{Float64}
    H = float.(hess)
    LinearAlgebra.ldiv!(invB, mcholesky!(H), 1)

    for i in 1:maxiter
        # выбор направления
        d = -invB * g
        
        # подбор шага вдоль d
        α = strong_backtracking(f, ∇f, x, d)
        if α == NaN
            return NewtonResult(false, x, i)  # α не найдено
        end
        
        # совершаем шаг и определяем s и y
        x_next .+= α*d
        g_next = ∇f(x_next)
        s = x_next - x
        y = g_next - g

        x = x_next

        # проверка сходимости
        g = g_next

        if LinearAlgebra.norm(g, 2) <= gtol
            return BFGSResult(true, x, 0)
        end

        # обновление обратного квази-гессиана для следующей итерации
        Bs = B*s
        invB .+= y*y' / dot(y, s) - Bs*Bs' / dot(s, Bs)
    end
    return BFGSResult(false, x, maxiter)
end