include("Points.jl")

using .Points

p = CartesianPoint(1, 1)
q = CartesianPoint(-2, 2)

println(p+q)
println(p-q)
println(norm(p))
println(arg(p)/pi)
println(Points.angle(p, q)/pi)
println(dot(p, q))
println(7 * p)
println(rotate(p, pi/4))

println("\n\n")

a = PolarPoint(p)
b = PolarPoint(-2*sqrt(2), -pi/4)

println(CartesianPoint(a+b))
println(CartesianPoint(a-b))
println(norm(a))
println(arg(a)/pi)
println(Points.angle(a, b)/pi)
println(dot(a, b))
println(CartesianPoint(7 * a))
println(CartesianPoint(rotate(a, pi/4)))

println("\n\n")

println(Points.angle(a, q)/pi)
println(dot(a, p))