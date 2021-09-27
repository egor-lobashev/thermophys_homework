include("pwlininterp.jl")

function eta_oc(M, sigma, e_over_k, T_min, T_max, filename)
    T_eta = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]
    f_eta = [1.0002, 1, 1, 1.0001, 1.0004, 1.0014, 1.0025, 1.0034, 1.0049, 1.0058, 1.0075]
    Omega_0 = [0.7 1.908 ; 0.75 1.841 ; 0.8 1.78 ; 0.85 1.725 ; 0.9 1.675 ; 0.95 1.627 ; 1 1.587 ; 1.05 1.549 ; 1.1 1.514 ; 1.15 1.482 ; 1.2 1.452 ; 1.25 1.424 ; 1.3 1.399 ; 1.35 1.375 ; 1.4 1.353 ; 1.45 1.333 ; 1.5 1.314 ; 1.55 1.296 ; 1.6 1.279 ; 1.65 1.264 ; 1.7 1.248 ; 1.75 1.234 ; 1.8 1.221 ; 1.85 1.209 ; 1.9 1.197 ; 1.95 1.186 ; 2 1.175 ; 2.1 1.156 ; 2.2 1.138 ; 2.3 1.122 ; 2.4 1.107 ; 2.5 1.093 ; 2.6 1.081 ; 2.7 1.069 ; 2.8 1.058 ; 2.9 1.048 ; 3 1.039 ; 3.1 1.03 ; 3.2 1.022 ; 3.3 1.014 ; 3.4 1.007 ; 3.5 0.9999 ; 3.6 0.9932 ; 3.7 0.987 ; 3.8 0.9811 ; 3.9 0.9755 ; 4 0.97 ; 4.1 0.9649 ; 4.2 0.96 ; 4.3 0.9553 ; 4.4 0.9507 ; 4.5 0.9464 ; 4.6 0.9422 ; 4.7 0.9382 ; 4.8 0.9343 ; 4.9 0.9305 ; 5 0.9269 ; 6 0.8963 ; 7 0.8727 ; 8 0.8538 ; 9 0.8379]

    Omega = pwlininterp(Omega_0[:, 1], Omega_0[:, 2])
    f_eta = pwlininterp(T_eta, f_eta)

    M /= 1000

    ys = T -> 8.44107*10^-5*sqrt(M * T * e_over_k) * f_eta(T) / (sigma^2*Omega(T))

    N = 100

    file = open(filename, "w")
    for j in 1:(N+1)
        T_1 = T_min + (T_max-T_min)/N * (j-1)
        # write(file, join([T_1, ", ", Omega(T_1/e_over_k), "\n"]))
        write(file, join([T_1, ", ", ys(T_1/e_over_k), "\n"]))
    end
end

eta_oc(44.009, 3.996, 190, 300, 1000, "CO_2.txt")
eta_oc(16.043, 3.822, 137, 100, 600, "CH_4.txt")
eta_oc(31.999, 3.433, 113, 100, 1000, "O_2.txt")