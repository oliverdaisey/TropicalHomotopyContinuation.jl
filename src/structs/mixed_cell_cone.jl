using Oscar

struct MixedCellConeFacet

    circuit::AbstractVector{Int}

end

struct MixedCellCone

    facets::AbstractVector{MixedCellConeFacet}

end

function facets(C::MixedCellCone)
    return C.facets
end

function Base.show(io::IO, C::MixedCellCone)

    print(io, "Mixed cell cone with $(length(facets(C))) facets")
end

function Base.show(io::IO, F::MixedCellConeFacet)

    print(io, "Mixed cell cone facet")
end

function mixed_cell_cone(facets::AbstractVector{MixedCellConeFacet})::MixedCellCone
    return MixedCellCone(facets)
end

"""
    mixed_cell_cone(candidate::MixedSupport, ambientSupport::MixedSupport)

Compute the mixed cell cone of `candidate` with ambient support `ambientSupport`.
"""
function mixed_cell_cone(candidate::MixedSupport, ambientSupport::MixedSupport)::MixedCellCone

    cayleyEmbedding = cayley_embedding(ambientSupport)
    
    for p in points(ambientSupport)
        if p in points(candidate)
            continue
        end
        # for all points p not in ambient support, get submatrix of cayleyEmbedding indexed by candidate and p
        
    end

end
