struct CayleyEmbedding

    supports::MixedSupport
    dim::Int

end

"""
    cayley_embedding(Δ::MixedSupport)::CayleyEmbedding

Construct the Cayley embedding of a mixed support `Δ`.
"""
function cayley_embedding(Δ::MixedSupport)::CayleyEmbedding

    dim = length(first(first(vector_of_points(Δ))))

    return CayleyEmbedding(Δ, dim)

end

function Base.show(io::IO, cayley::CayleyEmbedding)

    print(io, "Cayley embedding of $(length(supports(cayley))) point sets")
end

"""
    matrix(Δ::MixedSupport)

Given a mixed support `Δ`, returns its Cayley matrix.
"""
function matrix(Δ::MixedSupport)
    ptConfigurations = Vector{Vector{Int}}[vector_of_points(S) for S in supports(Δ)]
    numOfConfigurations = length(ptConfigurations)

    for (i, ptConfiguration) in enumerate(ptConfigurations)
        padding = zeros(Int, numOfConfigurations)
        padding[i] = 1
        ptConfigurations[i] = vcat.(ptConfiguration, Ref(padding))
    end

    return hcat(vcat(ptConfigurations...)...)
end

"""
    matrix(cayley::CayleyEmbedding)

Return the matrix representation of the Cayley embedding.
"""
function matrix(cayley::CayleyEmbedding)

    return matrix(supports(cayley))

end

"""
    supports(cayley::CayleyEmbedding)::MixedSupport

Return the mixed support involved in the Cayley embedding.
"""
function supports(cayley::CayleyEmbedding)::MixedSupport
    return cayley.supports
end

"""
    points(cayley::CayleyEmbedding)

Return the points of the Cayley embedding collected inside a vector.
"""
function points(cayley::CayleyEmbedding)
    return collect(eachcol(matrix(cayley)))
end

"""
    dim(cayley::CayleyEmbedding)

Return the dimension of the points defining the Cayley embedding. Note that this is not the same as the dimension of the Cayley embedding.
"""
function dim(cayley::CayleyEmbedding)
    return cayley.dim
end

"""
    Base.getindex(cayley::CayleyEmbedding, indices::MixedSupport)

Given a mixed support where each support is a subset of the supports of `cayley`, return the submatrix of `cayley` indexed by the supports in `indices`.
"""
function Base.getindex(cayley::CayleyEmbedding, indices::MixedSupport)

    @assert is_subset(indices, supports(cayley)) "The indices must be a subset of the supports of the Cayley embedding (got $indices, not a subset of $(supports(cayley)))"

    return matrix(indices)

end