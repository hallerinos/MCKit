function spinOrientation(ϕ::Float64,θ::Float64)
    x = sin(θ)cos(ϕ)
    y = sin(θ)sin(ϕ)
    z = cos(θ)
    return [x,y,z]
end