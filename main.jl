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
J = -0.5
𝐃₁ = [0,1,0]
𝐃₂ = [1,0,0]
B₀ = -0.3
𝐁 = [0,0,B₀]
# discretize the angles
N∠ = 100
Nϕ = 2*N∠
Nθ = N∠
# inverse temperature (larger values reject ΔE>0 with high probability)
β = 1e6
# how many Metropolis steps?
Nₛ = 1e5
# pbc or obc?
bc = "obc"

ps = generateParameters(Lx,Ly,J,𝐃₁,𝐃₂,𝐁,β,Nₛ,bc,Nϕ,Nθ,DMIHamiltonian)
configurations = initialConfigurations(ps)
# plot of the spin configurations
# p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
# display(plot(p...))

energy = []
convergence = β
iteration = 0
append!(energy,calculateEnergy(configurations, ps))
# perform the Metropolis algorithm as long as the energy fluctuation are larger than the average temperature fluctuations E≈β⁻¹
while convergence > β^(-1)
    configurations = metropolisAlgorithm(configurations, ps)
    append!(energy,calculateEnergy(configurations, ps))
    convergence = abs(energy[end]-energy[end-1])/(Lx*Ly)
    iteration += 1
    @info "perform Metropolis algorithm, iteration: " iteration convergence
end

##
p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
pE = plot(energy./(Lx*Ly),ylabel=L"\varepsilon_0",xlabel=L"N_s\times10^{%$(Int64(log10(Nₛ)))}",label=false,title="energy density")
pf = plot(p..., pE, layout=4)
savefig(pf,"classicalConfiguration.pdf")