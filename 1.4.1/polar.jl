function CartesianPoint(p::PolarPoint)
    return CartesianPoint(p.r*cos(p.phi), p.r*sin(p.phi))
end

function +(p::PolarPoint, q::PolarPoint)
    return PolarPoint(CartesianPoint(p) + CartesianPoint(q))
end

function -(p::PolarPoint, q::PolarPoint)
    return PolarPoint(CartesianPoint(p) - CartesianPoint(q))
end

function norm(p::PolarPoint)
    return p.r
end

function arg(p::PolarPoint)
    return p.phi
end

function *(alpha::Real, p::PolarPoint)
    return PolarPoint(alpha * p.r, p.phi)
end

function rotate(p::PolarPoint, phi::Real)
    return PolarPoint(p.r, (p.phi + phi)%(2*pi))
end