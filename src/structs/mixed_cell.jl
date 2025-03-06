@doc raw"""
    MixedCellCandidate

A candidate for a mixed cell, consisting of an active support and a chain of flats.
"""
struct MixedCell

    activeSupport::MixedSupport
    chainOfFlats::ChainOfFlats

end

@doc raw"""
    mixed_cell_candidate(activeSupport::MixedSupport, chainOfFlats::ChainOfFlats)::MixedCellCandidate

Construct a mixed cell candidate.
"""
function mixed_cell(activeSupport::MixedSupport, chainOfFlats::ChainOfFlats)::MixedCell
    return MixedCell(activeSupport, chainOfFlats)
end

@doc raw"""
    points(σ::MixedCellCandidate)

Return the points of the active support of the mixed cell candidate `σ`.
"""
function points(σ::MixedCell)
    return points(active_support(σ))
end

@doc raw"""
    chain_of_flats(σ::MixedCellCandidate)

Return the chain of flats of the mixed cell candidate `σ`.
"""
function chain_of_flats(σ::MixedCell)
    return σ.chainOfFlats
end

@doc raw"""
    active_support(σ::MixedCellCandidate)

Return the active support of the mixed cell candidate `σ`.
"""
function active_support(σ::MixedCell)
    return σ.activeSupport
end

@doc raw"""
    supports(σ::MixedCellCandidate)

Return the supports defining the active support of the mixed cell candidate `σ`.
"""
function supports(σ::MixedCell)
    return supports(active_support(σ))
end

function Base.length(σ::MixedCell)
    return length(supports(σ))
end

@doc raw"""
    Base.getindex(σ::MixedCellCandidate, p::Point)

Do not index mixed cell candidates directly as their data may be stale. Index the Tracker instead.
"""
function Base.getindex(σ::MixedCell, p::Point)
    error("Do not index mixed cells directly as their data may be stale. Index the Tracker instead.")
end

function Base.show(io::IO, σ::MixedCell)
    print(io, "Mixed cell with active points $(join(["$(p)" for p in points(σ)], ", ")) and chain of flats $(chain_of_flats(σ))")
end