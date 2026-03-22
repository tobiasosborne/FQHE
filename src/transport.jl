# transport.jl — Gaps → R_xx, R_xy vs B
#
# Semi-phenomenological transport model.
# Sources:
#   P10: Dykhne & Ruzin, PRB 50, 2369 (1994) — semicircle law
#   P11: Girvin lectures, §2-3 — general framework

export compute_Rxy, compute_Rxx

# TODO: Implement after gaps are computed (Phase C, Step 9)
function compute_Rxy(B_range, n_e, gaps, T; mat=GaAs)
    error("Not yet implemented — awaiting gap computations")
end

function compute_Rxx(B_range, n_e, gaps, T; mat=GaAs, μ_e=1e6*1e-4)
    error("Not yet implemented")
end
