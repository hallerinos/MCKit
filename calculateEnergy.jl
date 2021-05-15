function calculateEnergy(configurations, P)
    J = P["J"]
    𝐃₁ = P["𝐃₁"]
    𝐃₂ = P["𝐃₂"]
    bc = P["bc"]
    H = P["H"]
    𝐁 = P["𝐁"]
    
    Etot = 0
    # horizontal neighbor energies
    Etot += sum([H(configurations[x,y],configurations[x+1,y],J,𝐃₁,𝐁) for x=1:Lx-1 for y=1:Ly])
    # vertical neighbor energies
    Etot += sum([H(configurations[x,y],configurations[x,y+1],J,𝐃₂,𝐁) for x=1:Lx for y=1:Ly-1])
    if bc=="pbc"
        # horizontal boundary terms
        Etot += sum([H(configurations[Lx,y],configurations[1,y],J,𝐃₁,𝐁) for y=1:Ly])
        # vertical boundary terms
        Etot += sum([H(configurations[x,Ly],configurations[x,1],J,𝐃₂,𝐁) for x=1:Lx])
    end
    return Etot
end