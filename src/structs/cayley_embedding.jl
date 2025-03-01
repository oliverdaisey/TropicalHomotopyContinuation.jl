struct CayleyEmbedding

    supports::MixedSupport
    dim::Int
    matrix

end

function cayley_embedding(Δ::MixedSupport)::CayleyEmbedding
    ptConfigurations = Vector{Vector{Int}}[vector_of_points(S) for S in supports(Δ)]
    numOfConfigurations = length(ptConfigurations)

    for (i,ptConfiguration) in enumerate(ptConfigurations)
        padding = zeros(Int, numOfConfigurations)
        padding[i] = 1
        println("Before vcat: ", ptConfigurations[i])
        ptConfigurations[i] = vcat.(ptConfiguration, Ref(padding))
        println("After vcat: ", ptConfigurations[i])
    end


    dim = length(first(first(vector_of_points(Δ))))

    matrix = hcat(vcat(ptConfigurations...)...)
    return CayleyEmbedding(Δ, dim, matrix)
    
end

function Base.show(io::IO, cayley::CayleyEmbedding)

    print(io, "Cayley embedding of $(length(supports(cayley))) point sets")
end

function matrix(cayley::CayleyEmbedding)
    return cayley.matrix
end

function supports(cayley::CayleyEmbedding)
    return cayley.supports
end

function points(cayley::CayleyEmbedding)
    return collect(eachcol(matrix(cayley)))
end

function dim(cayley::CayleyEmbedding)
    return cayley.dim
end

function Base.getindex(cayley::CayleyEmbedding, point::Point, support::Support)

    i = findfirst(supports(cayley), support)
    M = matrix(cayleyEmbedding)
    N = dim(cayleyEmbedding)
    
    return i

end