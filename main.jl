using LinearAlgebra
using SparseArrays
using Plots, LaTeXStrings

include("calculateEnergy.jl")
include("spinOrientation.jl")
include("hamiltonian.jl")
include("metropolisAlgorithm.jl")
include("initialConfigurations.jl")
include("generateParameters.jl")

# system dimensions
size = 101
Lx = size
Ly = size
# coupling parameters
D = 1  # DMI amplitude
𝐃₁ = [0,D,0]  # horizontal DMI
𝐃₂ = [D,0,0]  # vertical DMI
J = -1.0*abs(D)  # exchange coupling
N∠ = 100  # discretize the angles
Nϕ = 2*N∠  # ϕ ranges from [0,2π]
Nθ = N∠  # θ ranges from [0,π]
β = 1e2  # inverse temperature
Nₛ = 1e6  # how many Metropolis steps per iteration?
bcs = ["obc"]  # boundary conditions
Nₘ = 1000  # maximum number of iterations
B₀s = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5].*(-abs(D))  # the magnetic fields which are simulated

##
for B₀ in B₀s, bc in bcs  # loop over the fields and boundary conditions
    fni = 0  # how many samples to simulate
    𝐁 = [0,0,B₀] 
    while fni<3
        ps = generateParameters(Lx,Ly,J,𝐃₁,𝐃₂,𝐁,β,Nₛ,bc,Nϕ,Nθ,DMIHamiltonian)  # simple Dict which contains everything
        configurations = initialConfigurations(ps)  # draw random config
        # plot of the spin configurations
        # p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
        # display(plot(p...))

        energy = []  # collect the configuration energies
        iteration = 1
        append!(energy,calculateEnergy(configurations, ps))  # append the initial energy
        # perform the Metropolis algorithm
        while Nₘ >= iteration
            # configurations, ΔEs, Nₐ = metropolisAlgorithm(configurations, ps)
            # E̅ = sum(ΔEs)/Nₛ
            metropolisAlgorithm!(configurations, ps)  # perform Nₛ MC steps
            E = calculateEnergy(configurations, ps)  # calculate energy of new config
            append!(energy,E)  # append energy
            Δε = (energy[end] - energy[end-1])/(Lx*Ly)  # compute the energy drift
            # nₐ = Nₐ/Nₛ*100  # step acceptance ratio
            @info "Running MC simulation for B=$B₀... iteration: " iteration E Δε
            iteration += 1
        end
        p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
        frac = Int64(floor(.9*iteration))  # show only a fraction of the energy history
        pE = plot(energy[end-frac:end]./(Lx*Ly),ylabel=L"\varepsilon_0",xlabel=L"N_s\times10^{%$(Int64(log10(Nₛ)))}",label=false,title="energy density")
        pf = plot(p..., pE, layout=4)
        # display(pf)

        path = "plots/B"*string(B₀)*"/"*"Lx"*string(Lx)
        if !isdir(path)
            mkpath("plots/B"*string(B₀)*"/")
        end
        sl = ["J",J,"D",D,bc,fni]
        fn = join([string(s) for s in sl])
        savefig(pf,path*"conf"*fn*".png")
        fni += 1
    end
end