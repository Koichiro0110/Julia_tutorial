# Conduct a basic DMRG calculation on the 1D XXZ Model
using ITensors
using HDF5 

let 
    # Define the number of sites in the chain
    # The random MPS must be defined in the same Hilbert space as the Hamiltonian
    number_of_sites = 100
    sites = siteinds("S=1/2", number_of_sites)
    ψ₀ = randomMPS(sites; linkdims = 2)              

    # Define the paramters used in the Hamiltonian
    J₁ = 1.0
    Δ = 1.5

    # Set up the Hamiltonian
    os = OpSum()

    # Set up the nearest-neighbor interactions
    for index = 1 : number_of_sites - 1
        os += 0.5 * J₁, "S+", index, "S-", index + 1
        os += 0.5 * J₁, "S-", index, "S+", index + 1
        os += Δ, "Sz", index, "Sz", index + 1
    end
    Hamiltonian = MPO(os, sites)

    # Parameters that control the DMRG optimization process
    nsweeps = 20
    maxdim = [20, 50, 200, 1000]
    cutoff = [1E-10]
    # Perform the DMRG calculation
    E₀, ψ = dmrg(Hamiltonian, ψ₀; nsweeps, maxdim, cutoff)
    
    println("The converged ground state energy is: $E₀")
    println("")
    # @show E₀

    # Measure one-point and two-point functions
    Sx = expect(ψ, "Sx"; sites = 1 : number_of_sites)
    Cxx = correlation_matrix(ψ, "Sx", "Sx"; sites = 1 : number_of_sites)
    @show Sx
    # @show Czz

    # Store the wavefunction and expectation value of observables into a HDF5 file
    h5open("Data/XXZ_Delta$(Δ)_N$(number_of_sites)_x.h5", "w") do file
        write(file, "E0", E₀)
        write(file, "Sx", Sx)
        write(file, "Cxx", Cxx)
        write(file, "Psi", ψ)
    end    

    return
end