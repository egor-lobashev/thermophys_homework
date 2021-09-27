module Points

export AbstractPoint, CartesianPoint, PolarPoint, norm, arg, rotate, angle, dot

import Base.:+, Base.:-, Base.:*

abstract type AbstractPoint end

struct CartesianPoint <: AbstractPoint
    x::Real
    y::Real
end

struct PolarPoint <: AbstractPoint
    r::Real
    phi::Real
    PolarPoint(r, phi) = r < 0 ? new(-r, (phi+pi) % (2*pi)) : new(r, phi % (2*pi))
end

include("cartesian.jl")
include("polar.jl")
include("abstract.jl")

end