function initialConfigurations(P,ctrl)
    dϕs = P["dϕs"]
    dθs = P["dθs"]
    Nϕ = P["Nϕ"]
    Nθ = P["Nθ"]

    # combine the angles to obtain the initial configuration
    if ctrl=="rand"
        # draw random angle configuration for each lattice position
        phis = [dϕs[rand(1:Nϕ)] for x=1:Lx, y=1:Ly]
        thetas = [dθs[rand(1:Nθ)] for x=1:Lx, y=1:Ly]
        configurations = [spinOrientation(phis[x,y],thetas[x,y]) for x=1:Lx, y=1:Ly]
    end
    if ctrl=="up"
        configurations = [spinOrientation(0.0,0.0) for x=1:Lx, y=1:Ly]
    end
    if ctrl=="down"
        configurations = [spinOrientation(0.0,1*π) for x=1:Lx, y=1:Ly]
    end

    return configurations
end