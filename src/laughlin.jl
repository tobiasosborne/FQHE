# laughlin.jl — Laughlin wavefunction overlaps (validation)
#
# Sources:
#   P2: Laughlin, PRL 50, 1395 (1983)
#   P8: Jain, Composite Fermions, Ch. 3

export laughlin_overlap

# TODO: Implement for validation (Phase C, Step 6 fallback)
# For first pass, validate via gap values instead.
function laughlin_overlap(ground_state, basis::FockBasis)
    error("Not yet implemented — validate via gaps first, implement if needed")
end
