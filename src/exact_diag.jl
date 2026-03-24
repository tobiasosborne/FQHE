# exact_diag.jl — Lanczos exact diagonalisation and gap extraction
#
# Sources:
#   KrylovKit.jl documentation
#   P04: Fano et al., PRB 34, 2670 (1986), Sec. IV

using KrylovKit, LinearAlgebra

export compute_spectrum, neutral_gap, charge_gap, GapTable, extrapolate_gap

"""
    compute_spectrum(H; nev=5) → (eigenvalues, eigenvectors)

Lowest `nev` eigenvalues/vectors of sparse Hermitian H.
Uses dense diag for dim ≤ 1500, Lanczos otherwise.
"""
function compute_spectrum(H; nev::Int=5)
    D = size(H, 1)
    nev = min(nev, D)

    if D <= 4000
        vals = eigvals(Matrix(H))
        return sort(real.(vals))[1:nev], nothing
    else
        vals, vecs, info = eigsolve(H, D > nev ? nev : D, :SR;
                                     issymmetric=true, tol=1e-12, maxiter=500)
        info.converged >= nev || @warn "Only $(info.converged)/$nev eigenvalues converged"
        return real.(vals[1:nev]), vecs[1:nev]
    end
end

"""
    neutral_gap(basis, VJ; nev=8) → (E₀, Δ_exc)

Ground state energy per particle and neutral excitation gap at fixed (N, twoS).
Δ_exc = E₁ − E₀ where E₁ is the lowest excited state.
Both in units of e²/(ε ℓ_B).
"""
function neutral_gap(basis::FockBasis, VJ::Vector{Float64}; nev::Int=8)
    H = build_hamiltonian(basis, VJ)
    vals, _ = compute_spectrum(H; nev=nev)
    E0 = vals[1]
    E0_per_particle = E0 / basis.N
    gap = length(vals) > 1 ? vals[2] - vals[1] : NaN
    return E0_per_particle, gap
end

"""
    charge_gap(N, ν) → (E₀/N, Δ_charge)

Charge gap Δ̃ = ε₊ + ε₋ via diagonalisation at 2S, 2S±1.
Handles parity: if Lz=0 sector is empty, uses the nearest available Lz.
"""
function charge_gap(N::Int, ν::Rational; nev::Int=5)
    twoS = sphere_flux(N, ν)
    E0 = _lowest_energy(N, twoS)

    # Quasihole: 2S + 1;  Quasiparticle: 2S − 1
    E_qh = _lowest_energy(N, twoS + 1)
    E_qp = _lowest_energy(N, twoS - 1)

    # Rescale to constant R = √(twoS/2)  (Fano convention)
    S = twoS / 2.0
    ε_minus = E_qh * sqrt((twoS + 1) / (2.0 * S)) - E0
    ε_plus  = E_qp * sqrt((twoS - 1) / (2.0 * S)) - E0
    return E0 / N, ε_plus + ε_minus
end

"""
    _lowest_energy(N, twoS) → E₀

Absolute ground state energy for N particles at monopole strength twoS.
Uses the smallest-|Lz| sector with correct parity.  By rotational invariance,
this sector contains components of ALL angular momentum multiplets L ≥ |Lz|,
so its lowest eigenvalue IS the absolute ground state energy.
"""
function _lowest_energy(N::Int, twoS::Int)
    # 2Lz must have same parity as N*twoS
    parity = mod(N * twoS, 2)
    twoLz_max = N * twoS - N * (N - 1)   # all N particles in highest orbitals
    VJ = coulomb_pseudopotentials(twoS)

    for twoLz in parity:2:twoLz_max
        basis = enumerate_fock_states(N, twoS; twoLz=twoLz)
        isempty(basis.states) && continue

        H = build_hamiltonian(basis, VJ)
        vals, _ = compute_spectrum(H; nev=1)
        return vals[1]
    end
    error("No states found for N=$N, twoS=$twoS")
end

"""Result of gap computation for one filling fraction."""
struct GapTable
    ν::Rational{Int}
    Ns::Vector{Int}
    E0_per_particle::Vector{Float64}
    gaps::Vector{Float64}          # charge gaps at each N
    gap_extrap::Float64            # extrapolated to N→∞
end

"""
    extrapolate_gap(Ns, gaps; order=0) → Δ(∞)

Finite-size extrapolation: Δ(N) = Δ(∞) + a₁/N + a₂/N² + ⋯
Least-squares fit of polynomial in 1/N.

`order` selects the highest power of 1/N:
  - `order=0`: auto-select (linear for ≤3 points, quadratic for ≥4).
  - `order=1`: linear 1/N only.
  - `order=2`: include 1/N² (captures even-odd shell effects).
"""
function extrapolate_gap(Ns::Vector{Int}, gaps::Vector{Float64}; order::Int=0)
    n = length(Ns)
    n >= 2 || return gaps[1]

    # Auto-select: quadratic when we have enough data, linear otherwise
    if order == 0
        order = n >= 4 ? 2 : 1
    end
    order = min(order, n - 1)   # can't fit more params than data points

    # Fit Δ(N) = c₀ + c₁/N + c₂/N² + ... via least squares
    x = 1.0 ./ Ns
    A = hcat([x .^ p for p in 0:order]...)
    c = A \ gaps
    return c[1]   # Δ(∞) = c₀
end
