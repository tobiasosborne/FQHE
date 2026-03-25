# dmrg_solver.jl — DMRG ground state and gap extraction on cylinder
#
# Wrapper around ITensors DMRG for FQHE at ν=1/3 on a cylinder.
# Initial state: CDW |100100100...⟩ (thin-cylinder ground state).
# Uses noise terms in early sweeps for mixing.
#
# Sources:
#   P25: Zaletel et al., PRL 110, 236801 (2013)

using ITensors, ITensorMPS

export dmrg_ground_state, dmrg_neutral_gap

"""
    dmrg_ground_state(Nphi, Ne, Ly; kwargs...) → (E0, psi, sites, H)

DMRG ground state for `Ne` electrons in `Nphi` LLL orbitals on a cylinder
of circumference `Ly`.

Returns ground state energy, MPS, site indices, and MPO (for reuse in gap calculation).
"""
function dmrg_ground_state(Nphi::Int, Ne::Int, Ly::Float64;
                            maxdim::Vector{Int}=[10, 20, 50, 100, 200, 400, 600, 800, 800, 800],
                            nsweeps::Int=length(maxdim),
                            cutoff::Float64=1e-10,
                            noise_scale::Float64=1e-3,
                            outputlevel::Int=0)
    # Build MPO
    sites, H, _ = build_cylinder_mpo(Nphi, Ne, Ly)

    # CDW initial state: |100100...⟩ (every 3rd orbital occupied for ν=1/3)
    state = [mod(i - 1, 3) == 0 ? "Occ" : "Emp" for i in 1:Nphi]
    # Verify particle count
    n_occ = count(s -> s == "Occ", state)
    if n_occ != Ne
        # Fallback: fill first Ne orbitals evenly spaced
        state = fill("Emp", Nphi)
        spacing = Nphi ÷ Ne
        for k in 0:Ne-1
            state[k * spacing + 1] = "Occ"
        end
    end
    psi0 = MPS(sites, state)

    # Noise schedule: start with noise, decay to zero
    noise = [noise_scale * 0.1^i for i in 0:min(4, nsweeps-1)]
    append!(noise, zeros(max(0, nsweeps - 5)))

    # Pad maxdim if shorter than nsweeps
    while length(maxdim) < nsweeps
        push!(maxdim, maxdim[end])
    end

    E0, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff=[cutoff], noise, outputlevel)

    return E0, psi, sites, H
end

"""
    dmrg_neutral_gap(Nphi, Ne, Ly; kwargs...) → (E0, E1, gap, psi0, psi1)

Compute neutral excitation gap via ground-state + excited-state DMRG.

The excited state is found by orthogonalizing against the ground state
with a penalty weight.
"""
function dmrg_neutral_gap(Nphi::Int, Ne::Int, Ly::Float64;
                           maxdim_gs::Vector{Int}=[10, 20, 50, 100, 200, 400, 600, 800, 800, 800],
                           maxdim_ex::Vector{Int}=[10, 20, 50, 100, 200, 400, 600, 800, 800, 800],
                           cutoff::Float64=1e-10,
                           noise_scale::Float64=1e-3,
                           weight::Float64=100.0,
                           outputlevel::Int=0)
    # Ground state
    E0, psi0, sites, H = dmrg_ground_state(Nphi, Ne, Ly;
                                             maxdim=maxdim_gs, cutoff, noise_scale, outputlevel)

    # Excited state: shifted CDW as initial guess
    state_ex = [mod(i, 3) == 0 ? "Occ" : "Emp" for i in 1:Nphi]
    n_occ = count(s -> s == "Occ", state_ex)
    if n_occ != Ne
        state_ex = fill("Emp", Nphi)
        spacing = Nphi ÷ Ne
        for k in 0:Ne-1
            state_ex[min(k * spacing + 2, Nphi)] = "Occ"
        end
    end
    psi_ex = MPS(sites, state_ex)

    nsweeps_ex = length(maxdim_ex)
    noise = [noise_scale * 0.1^i for i in 0:min(4, nsweeps_ex-1)]
    append!(noise, zeros(max(0, nsweeps_ex - 5)))
    while length(maxdim_ex) < nsweeps_ex
        push!(maxdim_ex, maxdim_ex[end])
    end

    E1, psi1 = dmrg(H, [psi0], psi_ex; nsweeps=nsweeps_ex, maxdim=maxdim_ex,
                      cutoff=[cutoff], noise, weight, outputlevel)

    gap = E1 - E0
    return E0, E1, gap, psi0, psi1
end
