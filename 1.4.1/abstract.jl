function angle(p::AbstractPoint, q::AbstractPoint)
    return abs(arg(p) - arg(q))
end

function dot(p::AbstractPoint, q::AbstractPoint)
    return norm(p) * norm(q) * cos(angle(p,q))
end