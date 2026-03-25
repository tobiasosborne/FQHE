# 05_itensors_smoke.jl — Smoke test: 4-site free-fermion DMRG vs exact diag
#
# H = -t Σ (c†_i c_{i+1} + h.c.) with t=1, N_sites=4, N_e=2 (half-filling)
# Exact spectrum: eigenvalues of 4×4 hopping matrix, fill lowest 2.

using ITensors, ITensorMPS
using LinearAlgebra

function main()
    println("=== ITensors Smoke Test: Free Fermions ===\n")

    N = 4
    Ne = 2
    t = 1.0

    # --- Exact diagonalization ---
    h = zeros(N, N)
    for i in 1:N-1
        h[i, i+1] = -t
        h[i+1, i] = -t
    end
    sp_energies = eigvals(Symmetric(h))
    E0_exact = sum(sp_energies[1:Ne])
    println("Single-particle energies: ", round.(sp_energies, digits=8))
    println("Exact GS energy (fill lowest $Ne): ", round(E0_exact, digits=10))

    # --- ITensors DMRG ---
    sites = siteinds("Fermion", N; conserve_qns=true)

    ampo = OpSum()
    for i in 1:N-1
        ampo += -t, "Cdag", i, "C", i+1
        ampo += -t, "Cdag", i+1, "C", i
    end
    H = MPO(ampo, sites)

    # Initial state: electrons on sites 1 and 3
    state = ["Occ", "Emp", "Occ", "Emp"]
    psi0 = MPS(sites, state)

    nsweeps = 5
    maxdim = [10, 20]
    cutoff = [1e-12]

    E0_dmrg, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, outputlevel=0)

    println("DMRG GS energy:                  ", round(E0_dmrg, digits=10))
    println()

    err = abs(E0_dmrg - E0_exact)
    println("Error: ", err)

    if err < 1e-10
        println("✓ SMOKE TEST PASSED (error < 1e-10)")
    else
        println("✗ SMOKE TEST FAILED")
        exit(1)
    end
end

main()
