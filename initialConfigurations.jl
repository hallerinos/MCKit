function initialConfigurations(P)
    dϕs = P["dϕs"]
    dθs = P["dθs"]
    Nϕ = P["Nϕ"]
    Nθ = P["Nθ"]

    # draw random angle configuration for each lattice position
    phis = [dϕs[rand(1:Nϕ)] for x=1:Lx, y=1:Ly]
    thetas = [dθs[rand(1:Nθ)] for x=1:Lx, y=1:Ly]
    # combine the angles to obtain the initial configuration
    configurations = [spinOrientation(phis[x,y],thetas[x,y]) for x=1:Lx, y=1:Ly]
end