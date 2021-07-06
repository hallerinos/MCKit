using LinearAlgebra
using SparseArrays
using Plots, LaTeXStrings
using JLD
using Logging

Logging.disable_logging(Logging.Info)  # use to suppress info messages

include("calculateEnergy.jl")
include("spinOrientation.jl")
include("hamiltonian.jl")
include("metropolisAlgorithm.jl")
include("initialConfigurations.jl")
include("generateParameters.jl")

# system dimensions
size = 4  # square lattice dimension
Lx = size  # x-dimension
Ly = size  # y-dimension
Nsamples = 100  # how many samples to simulate
# scale = 10
# Lx = Int64(floor(2560/scale))
# Ly = Int64(floor(1600/scale*sqrt(2)))
# coupling parameters
D = 1  # DMI amplitude
𝐃₁ = [0,D,0]  # horizontal DMI
𝐃₂ = [D,0,0]  # vertical DMI
J = -0.5*abs(D)  # exchange coupling
N∠ = 101  # discretize the angles
Nϕ = N∠  # ϕ ranges from [0,2π]
Nθ = N∠  # θ ranges from [0,π]
β = 1e3  # inverse temperature
Nₛ = floor(β*Lx*Ly/10)  # how many Metropolis steps per iteration?
bcs = ["obc"]  # boundary conditions
Nₘ = 24*5  # maximum number of iterations
ϵ = 0.01  # to which fraction the MC step may change the old local configuration
B₀s = (0:0.05:1.0).*(-abs(D))  # the magnetic fields which are simulated
save_history = false
save_config = false
do_plots = false

##
for B₀ in B₀s, bc in bcs  # loop over the fields and boundary conditions
    fni = 0  # how many samples to simulate
    𝐁 = [0,0,B₀] 
    while fni<Nsamples
        ps = generateParameters(Lx,Ly,J,𝐃₁,𝐃₂,𝐁,β,Nₛ,bc,Nϕ,Nθ,DMIHamiltonian)  # simple Dict which contains everything

        path = "results/B$B₀/Lx$Lx/"  # path to store stuff
        fn = "J$(J)D$(D)$(bc)$fni"  # string list of parameters
        if !isdir(path)
            mkpath(path)  # create the path if it does not exist
        end
        if isfile(path*fn*".jld")  # load configuration
            jlddata = load(path*fn*".jld")  # open config file
            configurations = jlddata["cfgs"]  # read config
        else
            configurations = initialConfigurations(ps,"rand")  # draw random config
        end
        # plot of the spin configurations
        # p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]
        # display(plot(p...))

        energy = []  # collect the configuration energies
        cfgs_hist = []
        iteration = 1
        append!(energy,calculateEnergy(configurations, ps)/(Lx*Ly))  # append the initial energy
        # perform the Metropolis algorithm
        if save_history
            open(path*fn*"_hist.dat", "w") do history_file
                write(history_file, "i x y S_x S_y S_z\n")
                [write(history_file, "$iteration $x $y ",string(configurations[x,y][1])," ",string(configurations[x,y][2])," ",string(configurations[x,y][3]), "\n") for x=1:Lx, y=1:Ly]
                        write(history_file, "\n")
            end
        end
        while Nₘ >= iteration
            # configurations, ΔEs, Nₐ = metropolisAlgorithm(configurations, ps)
            # E̅ = sum(ΔEs)/Nₛ
            metropolisAlgorithm!(configurations, ps, ϵ)  # perform Nₛ MC steps
            E = calculateEnergy(configurations, ps)  # calculate energy of new config
            ε = E/(Lx*Ly)
            append!(energy,ε)  # append energy
            ΔE = (energy[end] - energy[end-1])  # compute the energy drift
            # nₐ = Nₐ/Nₛ*100  # step acceptance ratio
            @info "Running MC simulation for B=$B₀... iteration: " iteration ε ΔE
            if save_history
                open(path*fn*"_hist.dat", "a") do history_file
                    [write(history_file, "$iteration $x $y ",string(configurations[x,y][1])," ",string(configurations[x,y][2])," ",string(configurations[x,y][3]), "\n") for x=1:Lx, y=1:Ly]
                    write(history_file, "\n")
                end
            end
            iteration += 1
            # p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]  # plot magnetization
            # # frac = Int64(floor(.9*iteration))  # show only a fraction of the energy history
            # pE = plot(energy,ylabel=L"\varepsilon_0",xlabel=L"N_s\times%$(Nₛ)",label=false,title="energy density")  # plot energy
            # pf = plot(p..., pE, layout=4)  # make a layout of 4 figures
            # savefig(pf,path*"conf"*fn*".png")  # figure
        end
        if save_config
            jldopen(path*fn*".jld", "w") do save_file
                write(save_file, "cfgs", configurations)
            end
        end

        open(path*fn*"_last.dat", "w") do history_file
            write(history_file, "i x y S_x S_y S_z\n")
            [write(history_file, "$iteration $x $y ",string(configurations[x,y][1])," ",string(configurations[x,y][2])," ",string(configurations[x,y][3]), "\n") for x=1:Lx, y=1:Ly]
        end

        if do_plots
            p = [heatmap(getindex.(configurations,i),title=L"S_%$i") for i=1:3]  # plot magnetization
            # frac = Int64(floor(.9*iteration))  # show only a fraction of the energy history
            pE = plot(energy,ylabel=L"\varepsilon_0",xlabel=L"N_s\times%$(Nₛ)",label=false,title="energy density")  # plot energy
            pf = plot(p..., pE, layout=4)  # make a layout of 4 figures
            savefig(pf,path*"conf"*fn*".png")  # figure
        end
        fni += 1
    end
end