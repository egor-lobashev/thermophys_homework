function +(p::CartesianPoint, q::CartesianPoint)
    return CartesianPoint(p.x + q.x, p.y + q.y)
end

function -(p::CartesianPoint, q::CartesianPoint)
    return CartesianPoint(p.x - q.x, p.y - q.y)
end

function norm(p::CartesianPoint)
    return sqrt(p.x^2 + p.y^2)
end

function arg(p::CartesianPoint)
    return p.y >= 0  ?  acos(p.x/norm(p))  :  2*pi - acos(p.x/norm(p))
end

function *(alpha::Real, p::CartesianPoint)
    return CartesianPoint(alpha * p.x, alpha * p.y)
end

function rotate(p::CartesianPoint, phi::Real)
    return CartesianPoint(cos(phi) * p.x - sin(phi) * p.y, sin(phi) * p.x + cos(phi) * p.y)
end

function PolarPoint(p::CartesianPoint)
    return PolarPoint(norm(p), arg(p))
end