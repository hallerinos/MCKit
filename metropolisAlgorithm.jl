function metropolisAlgorithm(configs,P)
    Nϕ = P["Nϕ"]
    Nθ = P["Nθ"]
    H = P["H"]
    β = P["β"]
    J = P["J"]
    𝐃₁ = P["𝐃₁"]
    𝐃₂ = P["𝐃₂"]
    dϕs = P["dϕs"]
    dθs = P["dθs"]

    bc = P["bc"]

    Lx = P["Lx"]
    Ly = P["Ly"]

    Nₛ = P["Nₛ"]

    cfgs = copy(configs)
    
    for s = 1:Nₛ
        # draw the position
        x = rand(1:Lx)
        y = rand(1:Ly)
        # draw the new angle
        ϕ = dϕs[rand(1:Nϕ)]
        θ = dθs[rand(1:Nθ)]
        # copy the old configuration
        oldS = cfgs[x,y]
        newS = spinOrientation(ϕ,θ)
        # compute the new and old local energies
        newε = 0
        oldε = 0
        if x<Lx
            newε += H(newS,cfgs[x+1,y],J,𝐃₁)
            oldε += H(oldS,cfgs[x+1,y],J,𝐃₁)
        end
        if x>1
            newε += H(cfgs[x-1,y],newS,J,𝐃₁)
            oldε += H(cfgs[x-1,y],oldS,J,𝐃₁)
        end
        if y<Ly
            newε += H(newS,cfgs[x,y+1],J,𝐃₂)
            oldε += H(oldS,cfgs[x,y+1],J,𝐃₂)
        end
        if y>1
            newε += H(cfgs[x,y-1],newS,J,𝐃₂)
            oldε += H(cfgs[x,y-1],oldS,J,𝐃₂)
        end
        # consider also boundary contributions
        if bc=="pbc"
            if x==1
                newε += H(cfgs[Lx,y],newS,J,𝐃₁)
                oldε += H(cfgs[Lx,y],oldS,J,𝐃₁)
            end
            if x==Lx
                newε += H(newS,cfgs[1,y],J,𝐃₁)
                oldε += H(oldS,cfgs[1,y],J,𝐃₁)
            end

            if y==1
                newε += H(cfgs[x,Ly],newS,J,𝐃₂)
                oldε += H(cfgs[x,Ly],oldS,J,𝐃₂)
            end
            if y==Lx
                newε += H(newS,cfgs[x,1],J,𝐃₂)
                oldε += H(oldS,cfgs[x,1],J,𝐃₂)
            end
        end
        # compute the energy difference
        ΔE = newε-oldε
        # acceptance probability
        pₐ = min(1, exp(-β*ΔE))
        # accept the new configuration if the acceptance probability pₐ is smaller than a random number
        if rand()>pₐ
            # @info "new configuration rejected"
            newS = oldS
            ΔE = 0
        end
        cfgs[x,y] = newS
    end

    return cfgs
end