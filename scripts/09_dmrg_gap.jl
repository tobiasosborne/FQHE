# 09_dmrg_gap.jl — Production DMRG gap for ν=1/3 on cylinder
#
# Finds the neutral excitation gap above the 3-fold degenerate ground state
# manifold (topological sectors of the Laughlin state).

using ITensors, ITensorMPS

push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using FQHE

include(joinpath(@__DIR__, "..", "src", "cylinder.jl"))
include(joinpath(@__DIR__, "..", "src", "dmrg_hamiltonian.jl"))

function make_cdw(sites, Nphi, Ne, shift)
    state = fill("Emp", Nphi)
    filled = 0
    for j in 0:Nphi-1
        if mod(j - shift, 3) == 0 && filled < Ne
            state[j+1] = "Occ"
            filled += 1
        end
    end
    for j in 1:Nphi
        filled >= Ne && break
        if state[j] == "Emp"
            state[j] = "Occ"
            filled += 1
        end
    end
    return MPS(sites, state)
end

function run_gap(Nphi, Ne, Ly, maxdim_schedule)
    Lx = cylinder_length(Nphi, Ly)
    println("  Nphi=$Nphi, Ne=$Ne, Ly=$Ly, Lx=$(round(Lx,digits=1)), Lx/Ly=$(round(Lx/Ly,digits=2))")

    sites, H, nt = build_cylinder_mpo(Nphi, Ne, Ly)
    println("  $nt terms, MPO bond dim=$(maxlinkdim(H))")

    nsweeps = length(maxdim_schedule)
    noise = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, fill(0.0, max(0, nsweeps - 5))...]
    noise = noise[1:nsweeps]

    # Best GS from 3 CDW shifts
    best_E = Inf
    best_psi = make_cdw(sites, Nphi, Ne, 0)
    for shift in 0:2
        psi0 = make_cdw(sites, Nphi, Ne, shift)
        E, psi = dmrg(H, psi0; nsweeps, maxdim=maxdim_schedule,
                       cutoff=[1e-10], noise, outputlevel=0)
        if E < best_E
            best_E = E
            best_psi = psi
        end
    end
    println("  E0 = $(round(best_E, digits=8)), E0/Ne = $(round(best_E/Ne, digits=6))")

    # Find states above GS manifold
    gs_list = MPS[best_psi]
    excitation_E = nothing

    for k in 1:4
        psi_init = make_cdw(sites, Nphi, Ne, k % 3)
        E, psi = dmrg(H, gs_list, psi_init; nsweeps, maxdim=maxdim_schedule,
                       cutoff=[1e-10], noise, weight=100.0, outputlevel=0)
        gap_k = E - best_E
        push!(gs_list, psi)

        if gap_k > 0.001  # significant gap — this is the excitation
            excitation_E = E
            println("  First excitation at state $k: E=$(round(E,digits=8)), gap=$(round(gap_k,digits=6))")
            break
        else
            println("  GS manifold state $k: E=$(round(E,digits=8)), δ=$(round(gap_k,digits=8))")
        end
    end

    if excitation_E === nothing
        println("  WARNING: no excitation found above GS manifold")
        return best_E, best_E, 0.0
    end

    gap = excitation_E - best_E
    return best_E, excitation_E, gap
end

function main()
    println("=== DMRG Neutral Gap for ν=1/3 ===\n")

    # Best configurations
    configs = [
        ("Ly=10, Nphi=21", 21, 7, 10.0,
         [10,20,50,100,200,400,600,800,800,800,800,800]),
        ("Ly=10, Nphi=24", 24, 8, 10.0,
         [10,20,50,100,200,400,600,800,800,800,800,800]),
        ("Ly=10, Nphi=30", 30, 10, 10.0,
         [10,20,50,100,200,400,800,1000,1000,1000,1000,1000]),
    ]

    results = []
    for (label, Nphi, Ne, Ly, maxdim) in configs
        println("\n--- $label ---")
        t0 = time()
        E0, E1, gap = run_gap(Nphi, Ne, Ly, maxdim)
        dt = time() - t0
        push!(results, (label=label, E0_N=E0/Ne, gap=gap, time=dt))
        println("  Time: $(round(dt, digits=1))s")
    end

    println("\n\n=== Summary ===")
    println("  Method              | E0/N      | Gap")
    println("  --------------------|-----------|--------")
    for r in results
        println("  $(rpad(r.label, 20))| $(lpad(round(r.E0_N,digits=6),9)) | $(round(r.gap,digits=6))")
    end
    println("  ED extrap (sphere)  |   0.085   | 0.0852")
    println("  Published DMRG      |   0.410   | 0.101")
end

main()
