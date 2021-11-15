module Integration

export Gauss, Kronrod, Midpoint, Trapezoid, Simpson
export integrate

abstract type AbstractMethod end

methods() = (
    Gauss,
    Kronrod,
    Midpoint,
    Trapezoid,
    Simpson,
)

"""
    integrate(f, a, b; method)

Интегрирование `f` на отрезке [`a`, `b`] методом `method`.

Список доступных методов: см. `Integration.methods()`
"""
integrate(f, a, b; method)

"""
    integrate(f, a, b, nnodes; method)

Интегрирование `f` на отрезке [`a`, `b`] составным вариантом метода `method`.
При этом отрезок предварительно разбивается `nnodes` равноотстоящими узлами.

Также см. [`integrate(f, a, b; method)`](@ref)
"""
integrate(f, a, b, nnodes; method)

# Диспетчеризация интерфейса
integrate(args...; method) = __integrate_impl(method, args...)

# Общий метод для составной формулы
function __integrate_impl(method::AbstractMethod, f, a, b, nnodes)
    x = range(a, b; length=nnodes-1)
    int = 0.0
    @views for (x₁, x₂) in zip(x[1:end-1], x[2:end])
        int += __integrate_impl(method, f, x₁, x₂)
    end
    return int
end


"Формула средних прямоугольников"
struct Midpoint <: AbstractMethod end

__integrate_impl(method::Midpoint, f, a, b) = (b-a) * f((b+a)/2)

# составной метод для формулы прямоугольников
function __integrate_impl(method::Midpoint, f, a, b, nnodes)
    h = (b - a) / (nnodes - 1)
    x = range(a + h/2, b; step=h)
    int = h * sum(f, x)
    return int
end

"Формула трапеций"
struct Trapezoid <: AbstractMethod end
__integrate_impl(method::Trapezoid, f, a, b) = (f(a) + f(b))/2 * (b-a)

"Формула Симпсона"
struct Simpson <: AbstractMethod end
__integrate_impl(method::Simpson, f, a, b) = (b-a)/6 * (f(a) + 4*f((a+b)/2) + f(b))
# __integrate_impl(method::Simpson, f, a, b, nnodes) = error("Не имплементирован")


# В качестве формулы Гаусса возьмите G_7
# https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula#Example
# Не забудьте пересчитать узлы и веса на отрезок [a, b]

function GK(f, a, b, nodes, weights)
    middle = (a+b)/2
    half = (b-a)/2
    ans = 0.0

    for i in 1:length(nodes)
        node = nodes[i]
        ans += f(middle + half*node) * weights[i] * half

        if node != 0
            ans += f(middle - half*node) * weights[i] * half
        end
    end

    return ans
end

struct Gauss <: AbstractMethod end
"Интгерирование методом Гаусса по 7 точкам."

G7_nodes = [0.949107912342759, 0.741531185599394, 0.405845151377397, 0]
G7_weights = [0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469]

__integrate_impl(method::Gauss, f, a, b) =
    GK(f, a, b, G7_nodes, G7_weights)

# В качестве формулы Кронрода возьмите K_15
# https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula#Example
# Не забудьте пересчитать узлы и веса на отрезок [a, b]

struct Kronrod <: AbstractMethod end

K15_nodes = [0.991455371120813, 0.949107912342759, 0.864864423359769, 0.741531185599394, 0.586087235467691, 0.405845151377397, 0.207784955007898, 0.0]
K15_weights = [0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267, 0.190350578064785, 0.204432940075298, 0.209482141084728]

__integrate_impl(method::Kronrod, f, a, b) =
    GK(f, a, b, K15_nodes,K15_weights)

end # module