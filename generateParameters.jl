function generateParameters(Lx,Ly,J,𝐃₁,𝐃₂,𝐁,β,Nₛ,bc,Nϕ,Nθ,H)
    P = Dict()

    P["Lx"] = Lx
    P["Ly"] = Ly
    
    P["J"] = J
    P["𝐃₁"] = 𝐃₁
    P["𝐃₂"] = 𝐃₂
    P["𝐁"] = 𝐁

    P["β"] = β
    P["Nₛ"] = Nₛ
    P["bc"] = bc

    P["Nϕ"] = Nϕ
    P["Nθ"] = Nθ
    δϕ = 2π/Nϕ
    δθ = π/Nθ
    dϕs = 0:δϕ:2π-δϕ
    dθs = 0:δθ:π-δθ
    P["dϕs"] = dϕs
    P["dθs"] = dθs

    P["H"] = H

    return P
end