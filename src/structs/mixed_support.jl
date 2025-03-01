using Oscar
export MixedSupport, vector_of_points, cayley_embedding, mixed_subdivision, supports

mutable struct MixedSupport

    supports::Tuple{Vararg{Support}} # collection of supports

end

function mixed_support(S::Tuple{Vararg{Support}})::MixedSupport
    return MixedSupport(S)
end

function supports(Δ::MixedSupport)::Tuple{Vararg{Support}}
    return Δ.supports
end

function Base.show(io::IO, Δ::MixedSupport)

    print(io, "Mixed support involving $(length(supports(Δ))) point supports")
end

# function cayley_embedding(Δ::MixedSupport)
#     ptConfigurations = Vector{Vector{Int}}[vector_of_points(S) for S in supports(Δ)]
#     numOfConfigurations = length(ptConfigurations)

#     for (i,ptConfiguration) in enumerate(ptConfigurations)
#         padding = zeros(Int, numOfConfigurations)
#         padding[i] = 1
#         ptConfigurations[i] = vcat.(ptConfiguration..., Ref(padding))
#     end

#     return hcat(vcat(ptConfigurations...)...)
    
# end

"""
    vector_of_points(Δ::MixedSupport)

Convenience function to return the points of a mixed support as a partitioned vector.
"""
function vector_of_points(Δ::MixedSupport)
    return [[p for p in points(s)] for s in supports(Δ)]
end

"""
    mixed_subdivision(Δ::MixedSupport)

Computes the mixed subdivision corresponding to the mixed support Δ.
"""
function mixed_subdivision(Δ::MixedSupport)
    minkowskiSum = minkowski_sum(points.(supports(Δ))...)
    mixedWeights = minkowski_sum(weights.(supports(Δ))...)

    subdivision = subdivision_of_points(convert.(Vector{Int}, minkowskiSum), convert.(QQFieldElem, mixedWeights))
    return subdivision
end