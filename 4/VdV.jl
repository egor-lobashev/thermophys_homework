include("all_roots.jl")

f(y) = (y*cosh(y) - sinh(y)) / (sinh(y)*cosh(y) - y)
g(y) = 1 + 2*f(y)*cosh(y) + f(y)^2

TX(y) = 27*f(y)*(f(y) + cosh(y)) / (4 * g(y)^2)
PX(y) = 27*f(y)^2*(1 - f(y)^2) / g(y)^2

Pr(V, T) = 8*T/(3*V - 1) - 3/V^2

function PX_VL_and_VG(T)
    TX_minus_T(y) = TX(y) - T
    y = all_roots(TX_minus_T, 0, 10, 1)[1]

    P = PX(y)
    Pr_minus_P(V) = Pr(V, T) - P

    VL_VG = all_roots(Pr_minus_P, 0.3334, 10, 3)
    return [P, VL_VG[1], VL_VG[3]]
end

for T in range(0.85, 1, step = 0.05)
    points = T==1 ? [1,1,1] : PX_VL_and_VG(T)
    println("T = ", T)
    println("PX = ", points[1], ", VL = ", points[2], ", VG = ", points[3], '\n')
    # println(points[2], ", ", points[1])
    # println(points[3], ", ", points[1])
end