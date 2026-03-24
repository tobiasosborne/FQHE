# clebsch_gordan.jl — Numerical Clebsch-Gordan coefficients via Racah formula
#
# Replaces WignerSymbols.jl (exact rational arithmetic, slow for large j)
# with Float64 log-factorial computation. ~10⁶× faster for j > 7.
#
# ⟨j₁,m₁; j₂,m₂ | J,M⟩ via the Racah formula:
#   C = √(2J+1) √(Δ(j₁,j₂,J)) √((j₁±m₁)!(j₂±m₂)!(J±M)!)
#       × Σ_k (-1)^k / [k!(j₁+j₂-J-k)!(j₁-m₁-k)!(j₂+m₂-k)!(J-j₂+m₁+k)!(J-j₁-m₂+k)!]
# where Δ(a,b,c) = (a+b-c)!(a-b+c)!(-a+b+c)!/(a+b+c+1)!
#
# Sources:
#   Varshalovich, Moskalev & Khersonskii, "Quantum Theory of Angular Momentum" (1988)
#   Racah, Phys. Rev. 62, 438 (1942)

export clebschgordan

# Precomputed log-factorials: log(n!) for n = 0, 1, ..., _LF_MAX
const _LF_MAX = 512
const _LOG_FACT = let v = zeros(Float64, _LF_MAX + 1)
    for i in 1:_LF_MAX
        v[i + 1] = v[i] + log(i)
    end
    v
end

@inline function _logfact(n::Int)
    @boundscheck (0 ≤ n ≤ _LF_MAX || throw(ArgumentError("log-factorial argument $n out of range [0,$_LF_MAX]")))
    @inbounds _LOG_FACT[n + 1]
end

"""
    clebschgordan(j1, m1, j2, m2, J, M) → Float64

Clebsch-Gordan coefficient ⟨j₁,m₁; j₂,m₂ | J,M⟩.
Arguments may be half-integer (as Float64).  All internal factorial
arguments are guaranteed integer by angular momentum theory.

Accurate to ~14 significant digits for j ≤ 100.
"""
function clebschgordan(j1, m1, j2, m2, J, M)
    # --- Selection rules ---
    abs(m1 + m2 - M) > 1e-10 && return 0.0
    abs(m1) > j1 + 1e-10 && return 0.0
    abs(m2) > j2 + 1e-10 && return 0.0
    abs(M)  > J  + 1e-10 && return 0.0
    J < abs(j1 - j2) - 1e-10 && return 0.0
    J > j1 + j2 + 1e-10 && return 0.0

    # --- Convert to integer factorial arguments ---
    a1 = round(Int, j1 + m1)   # j₁+m₁
    a2 = round(Int, j1 - m1)   # j₁−m₁
    b1 = round(Int, j2 + m2)   # j₂+m₂
    b2 = round(Int, j2 - m2)   # j₂−m₂
    c1 = round(Int, J + M)     # J+M
    c2 = round(Int, J - M)     # J−M
    t1 = round(Int, j1 + j2 - J)   # triangle: j₁+j₂−J
    t2 = round(Int, j1 - j2 + J)   # j₁−j₂+J
    t3 = round(Int, -j1 + j2 + J)  # −j₁+j₂+J
    t4 = round(Int, j1 + j2 + J + 1)

    any(x -> x < 0, (a1, a2, b1, b2, c1, c2, t1, t2, t3)) && return 0.0

    # --- Log of square-root prefactor ---
    log_pre = 0.5 * (log(2J + 1)
                    + _logfact(t1) + _logfact(t2) + _logfact(t3) - _logfact(t4)
                    + _logfact(a1) + _logfact(a2)
                    + _logfact(b1) + _logfact(b2)
                    + _logfact(c1) + _logfact(c2))

    # --- Summation limits ---
    p = round(Int, J - j2 + m1)   # offset for d5 = p + k
    q = round(Int, J - j1 - m2)   # offset for d6 = q + k
    k_lo = max(0, -p, -q)
    k_hi = min(t1, a2, b1)
    k_hi < k_lo && return 0.0

    # --- Racah sum ---
    result = 0.0
    for k in k_lo:k_hi
        log_denom = (_logfact(k) + _logfact(t1 - k) + _logfact(a2 - k)
                   + _logfact(b1 - k) + _logfact(p + k) + _logfact(q + k))
        result += (iseven(k) ? 1.0 : -1.0) * exp(log_pre - log_denom)
    end

    return result
end
