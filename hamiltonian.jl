# the Hamiltonian function of two local configurations S1 and S2
function DMIHamiltonian(S1::Vector, S2::Vector, J::Float64, 𝐃::Vector)
    ε = J*dot(S1,S2) + dot(𝐃, S1×S2) + dot(𝐁,S1+S2)
    return ε
end