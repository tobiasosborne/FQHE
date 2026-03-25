# 07_validate_mpo.jl — Validate ITensors MPO against sparse matrix ED
#
# For N_e=3, N_phi=9, Ly=12 (cylinder):
# Primary check: reconstruct full H matrix from MPO, compare with sparse ED
# Secondary check: DMRG energy (may not reach exact GS if degenerate)

using ITensors, ITensorMPS
using LinearAlgebra, Combinatorics

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using FQHE

include(joinpath(@__DIR__, "..", "src", "cylinder.jl"))
include(joinpath(@__DIR__, "..", "src", "dmrg_hamiltonian.jl"))

function build_sparse_cylinder_H(Ne::Int, Nphi::Int, Ly::Float64)
    states = collect(combinations(0:Nphi-1, Ne))
    D = length(states)
    state_to_idx = Dict(sort(s) => i for (i, s) in enumerate(states))
    H = zeros(Float64, D, D)
    for (β_idx, occ) in enumerate(states)
        for ii in 1:Ne-1, jj in ii+1:Ne
            a, b = occ[ii], occ[jj]
            H[β_idx, β_idx] += cylinder_matrix_element(a, b, a, b, Ly)
        end
        for ii in 1:Ne-1, jj in ii+1:Ne
            c, d = occ[ii], occ[jj]
            strip = setdiff(occ, [c, d])
            for a in 0:Nphi-1
                a in strip && continue
                for b in a+1:Nphi-1
                    b in strip && continue
                    a + b == c + d || continue
                    (a == c && b == d) && continue
                    V = cylinder_matrix_element(a, b, c, d, Ly)
                    abs(V) > 1e-15 || continue
                    new_occ = sort([strip; a; b])
                    α_idx = get(state_to_idx, new_occ, 0)
                    α_idx > 0 || continue
                    sgn = _fermion_sign(occ, c, d, a, b)
                    H[α_idx, β_idx] += V * sgn
                end
            end
        end
    end
    return Symmetric(H)
end

function _fermion_sign(occ, c, d, a, b)
    state = UInt64(0)
    for o in occ; state |= UInt64(1) << o; end
    sign = 1
    mask = (UInt64(1) << c) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    state &= ~(UInt64(1) << c)
    mask = (UInt64(1) << d) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    state &= ~(UInt64(1) << d)
    mask = (UInt64(1) << b) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    state |= UInt64(1) << b
    mask = (UInt64(1) << a) - 1
    isodd(count_ones(state & mask)) && (sign = -sign)
    return sign
end

function main()
    Ne = 3
    Nphi = 9
    Ly = 12.0

    println("=== Phase 2: MPO vs Sparse Matrix Validation ===")
    println("  N_e=$Ne, N_phi=$Nphi, Ly=$Ly\n")

    # --- (a) Sparse ED ---
    println("Building sparse Hamiltonian...")
    H_sparse = build_sparse_cylinder_H(Ne, Nphi, Ly)
    vals_ed = eigvals(Matrix(H_sparse))
    println("  Hilbert dim: $(size(H_sparse, 1))")
    println("  ED spectrum (first 6): ", round.(vals_ed[1:6], digits=8))

    # --- (b) Build MPO and reconstruct full matrix ---
    println("\nBuilding MPO...")
    sites, H_mpo, n_terms = build_cylinder_mpo(Nphi, Ne, Ly)
    println("  OpSum terms: $n_terms")

    println("  Reconstructing full matrix from MPO...")
    all_states = collect(combinations(0:Nphi-1, Ne))
    D = length(all_states)
    mps_list = Vector{MPS}(undef, D)
    for (idx, occ) in enumerate(all_states)
        st = fill("Emp", Nphi)
        for o in occ; st[o+1] = "Occ"; end
        mps_list[idx] = MPS(sites, st)
    end

    H_from_mpo = zeros(D, D)
    for i in 1:D, j in i:D
        H_from_mpo[i,j] = inner(mps_list[i]', H_mpo, mps_list[j])
        H_from_mpo[j,i] = H_from_mpo[i,j]
    end

    vals_mpo = eigvals(Symmetric(H_from_mpo))
    println("  MPO spectrum (first 6): ", round.(vals_mpo[1:6], digits=8))

    # --- Primary check: full matrix comparison ---
    max_diff = maximum(abs.(H_from_mpo - Matrix(H_sparse)))
    println("\n  Max |H_MPO - H_ED| = $max_diff")

    if max_diff < 1e-10
        println("  ✓ FULL MATRIX VALIDATION PASSED (max diff < 1e-10)")
    else
        println("  ✗ FULL MATRIX VALIDATION FAILED")
        exit(1)
    end

    # --- Secondary check: DMRG (informational) ---
    println("\nRunning DMRG (informational, degenerate GS expected)...")
    state = [mod(i - 1, 3) == 0 ? "Occ" : "Emp" for i in 1:Nphi]
    psi0 = MPS(sites, state)

    nsweeps = 20
    maxdim = fill(200, nsweeps)
    cutoff = [1e-12]
    noise = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, fill(0.0, nsweeps - 5)...]

    E0_dmrg, psi = dmrg(H_mpo, psi0; nsweeps, maxdim, cutoff, noise, outputlevel=0)

    println("  E0 (DMRG):   $E0_dmrg")
    println("  E0 (ED):     $(vals_ed[1])")
    println("  E2 (ED):     $(vals_ed[3])  (next level above GS doublet)")
    err = abs(E0_dmrg - vals_ed[1])
    println("  |E0_DMRG - E0_ED| = $err")

    if err < 1e-6
        println("  ✓ DMRG found exact ground state")
    elseif abs(E0_dmrg - vals_ed[3]) < 1e-6
        println("  Note: DMRG found E2 (GS is doubly degenerate, DMRG converged to next singlet)")
        println("  This is a known issue with degenerate CDW states on small cylinders")
    else
        println("  Note: DMRG converged to E = $E0_dmrg (may need more sweeps/noise)")
    end

    println("\n=== MPO construction VALIDATED ===")
end

main()
