@doc raw"""
    MixedCellCandidate

A candidate for a mixed cell, consisting of an active support and a chain of flats.
"""
struct MixedCell

    activeSupport::MixedSupport
    chainOfFlats::ChainOfFlats

end

@doc raw"""
    mixed_cell(activeSupport::MixedSupport, chainOfFlats::ChainOfFlats)::MixedCellCandidate

Construct a mixed cell candidate.
"""
function mixed_cell(activeSupport::MixedSupport, chainOfFlats::ChainOfFlats)::MixedCell
    return MixedCell(activeSupport, chainOfFlats)
end

@doc raw"""
    points(σ::MixedCell)

Return the points of the active support of the mixed cell candidate `σ`.
"""
function points(σ::MixedCell)
    return points(active_support(σ))
end

@doc raw"""
    chain_of_flats(σ::MixedCell)

Return the chain of flats of the mixed cell candidate `σ`.
"""
function chain_of_flats(σ::MixedCell)
    return σ.chainOfFlats
end

@doc raw"""
    active_support(σ::MixedCell)

Return the active support of the mixed cell candidate `σ`.
"""
function active_support(σ::MixedCell)
    return σ.activeSupport
end

@doc raw"""
    supports(σ::MixedCell)

Return the supports defining the active support of the mixed cell candidate `σ`.
"""
function supports(σ::MixedCell)
    return supports(active_support(σ))
end

function Base.length(σ::MixedCell)
    return length(supports(σ))
end

@doc raw"""
    Base.getindex(σ::MixedCell, p::Point)

Do not index mixed cell candidates directly as their data may be stale. Index the Tracker instead.
"""
function Base.getindex(σ::MixedCell, p::Point)
    error("Do not index mixed cells directly as their data may be stale. Index the Tracker instead.")
end

function Base.show(io::IO, σ::MixedCell)
    print(io, "Mixed cell with active points $(join(["$(p)" for p in points(σ)], ", ")) and chain of flats $(chain_of_flats(σ))")
end

@doc raw"""
    transform_linear_support(σ::MixedCell)

Returns the mixed support obtained by replacing the chain of flats of `σ` with a tuple of hypersurfaces whose intersection equals the span of the corresponding fine structure cone.`
"""
function transform_linear_support(σ::MixedCell)::MixedSupport
    
    Δ = active_support(σ)
    gens = indicator_vector.(flats(chain_of_flats(σ)))
    push!(gens, ones(Int, length(first(gens))))
    # find basis of orthogonal complement of span(gens)
    gensMatrix = Oscar.matrix(gens)
    _, kernel = Oscar.nullspace(gensMatrix)
    # convert kernel to list of vectors
    kernel = [kernel[:, i] for i in 1:ncols(kernel)]
    # add the vectors to pts
    pts = [point(g) for g in kernel]

    S = Support[]
    for p in pts
        push!(S, support([point(zeros(Int, length(p))), p], [0, 0]))
    end

    println(pts)
    return mixed_support((S..., supports(σ)...))
end