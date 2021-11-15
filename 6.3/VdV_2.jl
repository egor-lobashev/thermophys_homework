function norm(x)
    return sqrt(x[1]^2 + x[2]^2)
end

function newtonsys(f, x, J; maxiter=50, xtol=1e-6, ftol=1e-6)
    x = float(copy(x))
    δx, y = similar.((x, x))
    for i in 1:maxiter
        y .= f(x)
        δx .= .- (J(x) \ y)
        x .+= δx
        if norm(δx) < xtol || norm(y) < ftol
            return x
        end
    end
    error("Превышено число итераций.")
end

function P_mu(V, T)
    return [
        8*T/(3*V - 1) - 3/V^2,
        -T*log((3*V - 1)/(2*exp(-0.5))) + T/(3*V - 1) - 9/(4*V),
    ]
end

function P_mu_derivative(V, T)
    return [
        -24*T/(3*V-1)^2+6/V^3;
        -T/(V-1/3)-3*T/(3*V-1)^2+9/(4*V^2)
    ]
end

function f(x, T)
    VL, VG = x
    return P_mu(VL, T) - P_mu(VG, T)
end

function J(x, T)
    VL, VG = x
    return [P_mu_derivative(VL, T)      -P_mu_derivative(VG, T)]
end

for T in range(0.85, 0.95, step = 0.02)
    VL_VG = newtonsys(x->f(x, T), [0.5, 4], x->J(x, T))
    # println("T = ", T)
    # print("P = ", P_mu(VL_VG[1], T)[1])
    # println(", VL = ", VL_VG[1], ", VG = ", VL_VG[2], '\n')
    println(VL_VG[1], ", ", P_mu(VL_VG[1], T)[1])
    println(VL_VG[2], ", ", P_mu(VL_VG[2], T)[1])
end