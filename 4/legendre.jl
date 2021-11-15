include("all_roots.jl")

P_7 = x -> 1/16 * (427*x^7 - 693*x^5 + 315x^3 - 35*x)
roots = all_roots(P_7, -1, 1, 7)

for r in roots
    println(r)
end