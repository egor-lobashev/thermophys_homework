include("cholesky.jl")
include("backtracking.jl")

struct NewtonResult{T<:Real}
    converged::Bool      # метод сошёлся
    argument::Vector{T}  # найденный минимум
    iters::Int           # число прошедших итераций
end

function newton(
    f::Function,
    ∇f::Function,
    hess::Function,
    x0::AbstractVector;
    gtol::Real=1e-5,
    maxiter::Integer=200,
)
    x = float.(x0)  # копирование и приведение к флоат-числам

    # проверяем сходимость по норме-2 градиента
    g = ∇f(x)

    if LinearAlgebra.norm(g, 2) <= gtol
        return NewtonResult(true, x, 0)
    end

    d = Vector{Float64}(undef, length(x))

    for i in 1:maxiter
        # выбор направления d из модифицированного Гессиана Hᵐ
        Hᵐ = mcholesky!(float.(hess(x)))
        LinearAlgebra.ldiv!(d, Hᵐ, -g)

        # выбор шага вдоль d
        α = strong_backtracking(f, ∇f, x, d)
        if α == NaN
            return NewtonResult(false, x, i)  # α не найдено
        end
        
        # совершаем шаг
        x .+= α*d

        # проверяем сходимость
        g = ∇f(x)

        if LinearAlgebra.norm(g, 2) <= gtol
            return NewtonResult(true, x, i)
        end
    end
    return NewtonResult(false, x, maxiter)
end