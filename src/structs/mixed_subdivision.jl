using Oscar
export MixedSubdivision, vector_of_points, cayley_embedding

mutable struct MixedSubdivision

    supports::Tuple{Vararg{Support}} # collection of supports

end

function mixed_subdivision(S::Tuple{Vararg{Support}})::MixedSubdivision
    return MixedSubdivision(S)
end

function supports(Δ::MixedSubdivision)::Tuple{Vararg{Support}}
    return Δ.supports
end

function Base.show(io::IO, Δ::MixedSubdivision)

    print(io, "Mixed subdivision of $(length(supports(Δ))) point supports")
end

function cayley_embedding(Δ::MixedSubdivision)
    pts = [vector_of_points(S) for S in supports(Δ)]

    return pts
    
end

"""
    vector_of_points(Δ::MixedSubdivision)

Convenience function to return the points of a mixed subdivision as a vector.
"""
function vector_of_points(Δ::MixedSubdivision)
    return [p for s in supports(Δ) for p in points(s)]
end