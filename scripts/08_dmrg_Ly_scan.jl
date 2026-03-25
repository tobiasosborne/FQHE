# 08_dmrg_Ly_scan.jl — DMRG gap vs cylinder circumference
#
# Strategy: the s=0 (q_y=0) diagonal terms give a large constant shift that depends
# on Ne but NOT on the state. We compute the Madelung shift exactly for the CDW
# configuration and subtract it, so DMRG works with a well-scaled Hamiltonian.
# The gap is unaffected.

using ITensors, ITensorMPS
using LinearAlgebra

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using FQHE

include(joinpath(@__DIR__, "..", "src", "cylinder.jl"))
include(joinpath(@__DIR__, "..", "src", "dmrg_hamiltonian.jl"))
include(joinpath(@__DIR__, "..", "src", "dmrg_solver.jl"))

function make_cdw_state(sites, Nphi, Ne, shift)
    """CDW product state with given shift: |0^shift 1 0 0 1 0 0 ...⟩"""
    state = fill("Emp", Nphi)
    filled = 0
    for j in 0:Nphi-1
        if mod(j - shift, 3) == 0 && filled < Ne
            state[j+1] = "Occ"
            filled += 1
        end
    end
    # Fallback: fill remaining
    for j in 1:Nphi
        filled >= Ne && break
        if state[j] == "Emp"
            state[j] = "Occ"
            filled += 1
        end
    end
    return MPS(sites, state)
end

function run_dmrg_gap(Nphi::Int, Ne::Int, Ly::Float64, maxdim_vals::Vector{Int})
    println("  Building MPO (Nphi=$Nphi, Ne=$Ne, Ly=$Ly)...")
    sites, H, n_terms = build_cylinder_mpo(Nphi, Ne, Ly)
    println("    $n_terms OpSum terms")

    nsweeps = length(maxdim_vals)
    noise = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, fill(0.0, max(0, nsweeps - 5))...]
    noise = noise[1:nsweeps]

    # At ν=1/3, the GS manifold is 3-fold degenerate (topological sectors).
    # Find all 3 ground states, then the first excitation above the manifold.
    gs_list = MPS[]
    energies = Float64[]

    for k in 0:3  # 3 ground states + 1 excitation
        psi_init = make_cdw_state(sites, Nphi, Ne, k)
        if isempty(gs_list)
            E, psi = dmrg(H, psi_init; nsweeps, maxdim=maxdim_vals,
                           cutoff=[1e-10], noise, outputlevel=0)
        else
            E, psi = dmrg(H, gs_list, psi_init; nsweeps, maxdim=maxdim_vals,
                           cutoff=[1e-10], noise, weight=100.0, outputlevel=0)
        end
        push!(gs_list, psi)
        push!(energies, E)
        println("    State $k: E = $(round(E, digits=8))")
    end

    E0 = energies[1]
    gap = energies[4] - E0
    println("    Gap (E4-E1) = $gap")
    return E0, energies[4], gap
end

function main()
    println("=== DMRG Gap vs Cylinder Circumference ===\n")

    # For each config, Nphi/3 electrons, aspect ratio ~1.5-2
    configs = [
        (8.0,  15, [10,20,50,100,200,200,200,200,200,200,200,200]),
        (10.0, 21, [10,20,50,100,200,400,400,400,400,400,400,400]),
        (12.0, 27, [10,20,50,100,200,400,800,800,800,800,800,800]),
    ]

    println("  Ly   | Nphi | Ne | Lx    | Lx/Ly | E0/Ne      | gap        | time(s)")
    println("  -----|------|----| ------|-------|------------|------------|-------")

    for (Ly, Nphi, maxdim) in configs
        Ne = Nphi ÷ 3
        Lx = cylinder_length(Nphi, Ly)
        aspect = Lx / Ly

        t0 = time()
        E0, E1, gap = run_dmrg_gap(Nphi, Ne, Ly, maxdim)
        dt = time() - t0

        println("  $(lpad(Ly, 4)) | $(lpad(Nphi, 4)) | $(lpad(Ne, 2)) | $(lpad(round(Lx,digits=1),5)) | $(lpad(round(aspect,digits=2),5)) | $(lpad(round(E0/Ne, digits=6),10)) | $(lpad(round(gap, digits=6),10)) | $(round(dt, digits=1))")
    end

    println("\n  Ref: E0/N = -0.4102 (thermo), gap = 0.101 (DMRG lit), 0.085 (ED extrap)")
end

main()
