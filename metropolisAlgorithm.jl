function metropolisAlgorithm(configs,P)
    Nϕ = Int64(P["Nϕ"])
    Nθ = Int64(P["Nθ"])
    H = P["H"]
    β = P["β"]
    J = P["J"]
    𝐁 = P["𝐁"]
    𝐃₁ = P["𝐃₁"]
    𝐃₂ = P["𝐃₂"]
    dϕs = P["dϕs"]
    dθs = P["dθs"]

    bc = P["bc"]

    Lx = Int64(P["Lx"])
    Ly = Int64(P["Ly"])

    Nₛ = Int64(P["Nₛ"])

    cfgs = copy(configs)
    
    ΔEs = Vector{Float64}(undef, Nₛ)
    Nₐ = 0
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
            newε += H(newS,cfgs[x+1,y],J,𝐃₁,𝐁)
            oldε += H(oldS,cfgs[x+1,y],J,𝐃₁,𝐁)
        end
        if x>1
            newε += H(cfgs[x-1,y],newS,J,𝐃₁,𝐁)
            oldε += H(cfgs[x-1,y],oldS,J,𝐃₁,𝐁)
        end
        if y<Ly
            newε += H(newS,cfgs[x,y+1],J,𝐃₂,𝐁)
            oldε += H(oldS,cfgs[x,y+1],J,𝐃₂,𝐁)
        end
        if y>1
            newε += H(cfgs[x,y-1],newS,J,𝐃₂,𝐁)
            oldε += H(cfgs[x,y-1],oldS,J,𝐃₂,𝐁)
        end
        # consider also boundary contributions
        if bc=="pbc"
            if x==1
                newε += H(cfgs[Lx,y],newS,J,𝐃₁,𝐁)
                oldε += H(cfgs[Lx,y],oldS,J,𝐃₁,𝐁)
            end
            if x==Lx
                newε += H(newS,cfgs[1,y],J,𝐃₁,𝐁)
                oldε += H(oldS,cfgs[1,y],J,𝐃₁,𝐁)
            end

            if y==1
                newε += H(cfgs[x,Ly],newS,J,𝐃₂,𝐁)
                oldε += H(cfgs[x,Ly],oldS,J,𝐃₂,𝐁)
            end
            if y==Lx
                newε += H(newS,cfgs[x,1],J,𝐃₂,𝐁)
                oldε += H(oldS,cfgs[x,1],J,𝐃₂,𝐁)
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
            Nₐ -= 1
        end
        ΔEs[s] = ΔE
        cfgs[x,y] = newS
        Nₐ += 1
    end

    return cfgs, ΔEs, Nₐ
end

function metropolisAlgorithm!(cfgs,P)
    Nϕ = Int64(P["Nϕ"])
    Nθ = Int64(P["Nθ"])
    H = P["H"]
    β = P["β"]
    J = P["J"]
    𝐁 = P["𝐁"]
    𝐃₁ = P["𝐃₁"]
    𝐃₂ = P["𝐃₂"]
    dϕs = P["dϕs"]
    dθs = P["dθs"]

    bc = P["bc"]

    Lx = Int64(P["Lx"])
    Ly = Int64(P["Ly"])

    Nₛ = Int64(P["Nₛ"])

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
            newε += H(newS,cfgs[x+1,y],J,𝐃₁,𝐁)
            oldε += H(oldS,cfgs[x+1,y],J,𝐃₁,𝐁)
        end
        if x>1
            newε += H(cfgs[x-1,y],newS,J,𝐃₁,𝐁)
            oldε += H(cfgs[x-1,y],oldS,J,𝐃₁,𝐁)
        end
        if y<Ly
            newε += H(newS,cfgs[x,y+1],J,𝐃₂,𝐁)
            oldε += H(oldS,cfgs[x,y+1],J,𝐃₂,𝐁)
        end
        if y>1
            newε += H(cfgs[x,y-1],newS,J,𝐃₂,𝐁)
            oldε += H(cfgs[x,y-1],oldS,J,𝐃₂,𝐁)
        end
        # consider also boundary contributions
        if bc=="pbc"
            if x==1
                newε += H(cfgs[Lx,y],newS,J,𝐃₁,𝐁)
                oldε += H(cfgs[Lx,y],oldS,J,𝐃₁,𝐁)
            end
            if x==Lx
                newε += H(newS,cfgs[1,y],J,𝐃₁,𝐁)
                oldε += H(oldS,cfgs[1,y],J,𝐃₁,𝐁)
            end

            if y==1
                newε += H(cfgs[x,Ly],newS,J,𝐃₂,𝐁)
                oldε += H(cfgs[x,Ly],oldS,J,𝐃₂,𝐁)
            end
            if y==Lx
                newε += H(newS,cfgs[x,1],J,𝐃₂,𝐁)
                oldε += H(oldS,cfgs[x,1],J,𝐃₂,𝐁)
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
        end
        cfgs[x,y] = newS
    end
end