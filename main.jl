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
size = 11
Lx = size
Ly = size
# coupling parameters
D = 1.0
𝐃₁ = [0,D,0]
𝐃₂ = [D,0,0]
J = -0.5*abs(D)
# discretize the angles
N∠ = 10000
Nϕ = 2*N∠  # ϕ ranges from [0,2π]
Nθ = N∠  # θ ranges from [0,π]
# inverse temperature (larger values reject ΔE>0 with high probability)
β = 1e2
# how many Metropolis steps per iteration?
Nₛ = 1e5
bc = "obc"

Nₘ = 100  # maximum number of iterations

B₀s = [0.0,0.1, 0.2, 0.3, 0.4].*(-abs(D))  # the magnetic fields which are simulated

for B₀ in B₀s, bc in ["obc"]  # loop over the fields and boundary conditions
    fni = 0  # how many samples to simulate
    𝐁 = [0,0,B₀] 
    while fni<3
        ps = generateParameters(Lx,Ly,J,𝐃₁,𝐃₂,𝐁,β,Nₛ,bc,Nϕ,Nθ,DMIHamiltonian)
        configurations = initialConfigurations(ps)
        # plot of the spin configurations
        # p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
        # display(plot(p...))

        energy = []
        iteration = 1
        Nₐ = 2
        append!(energy,calculateEnergy(configurations, ps))
        # perform the Metropolis algorithm
        while Nₘ >= iteration && Nₐ > 1
            configurations, ΔEs, Nₐ = metropolisAlgorithm(configurations, ps)
            E̅ = sum(ΔEs)/Nₛ
            append!(energy,calculateEnergy(configurations, ps))
            Δε = (energy[end] - energy[end-1])/(Lx*Ly)
            nₐ = Nₐ/Nₛ*100
            @info "Running MC simulation for B=$B₀... iteration: " iteration E̅ Δε nₐ
            iteration += 1
        end
        p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
        frac = Int64(floor(.5*iteration))
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