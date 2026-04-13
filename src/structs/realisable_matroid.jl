@doc raw"""
    RealisableMatroid

An interface to a realisable matroid that provides access to the realisation matrix.
"""
struct RealisableMatroid
    realisationMatrix::MatElem{<:FieldElem}
    rank::Int
    zeroColumns::Set{Int}
    nonzeroColumns::Vector{Int}
    rankCache::Dict{Vector{Int}, Int}
end

@doc raw"""
    matrix(M::RealisableMatroid)

Return the realisation matrix of the realisable matroid `M`.
"""
function matrix(M::RealisableMatroid)
    return M.realisationMatrix
end

@doc raw"""
    matroid(M::RealisableMatroid)

Return the underlying Oscar matroid of the realisable matroid `M`.
"""
function matroid(M::RealisableMatroid)
    return matroid_from_matrix_columns(matrix(M))
end

@doc raw"""
    base_field(M::RealisableMatroid)

Return the field over which `M` is realised.
"""
function base_field(M::RealisableMatroid)
    return base_ring(matrix(M))
end

function Base.show(io::IO, M::RealisableMatroid)
    print(io, "Realisable matroid of rank ", rank(M), " with ground set ", ground_set(M))
end

@doc raw"""
    zero_columns(M::RealisableMatroid)

Return the set of column indices of `M` that are entirely zero.
"""
zero_columns(M::RealisableMatroid) = M.zeroColumns

@doc raw"""
    nonzero_columns(M::RealisableMatroid)

Return the sorted vector of column indices of `M` that are not entirely zero.
"""
nonzero_columns(M::RealisableMatroid) = M.nonzeroColumns

@doc raw"""
    matroid(A::MatElem{<:FieldElem})

Construct a realisable matroid from the realisation matrix `A`.
"""
function matroid(A::MatElem{<:FieldElem})
    r = Oscar.rank(A)
    n = ncols(A)
    zero_cols = Set{Int}()
    nonzero_cols = Int[]
    for j in 1:n
        is_zero_col = true
        for i in 1:nrows(A)
            if !iszero(A[i, j])
                is_zero_col = false
                break
            end
        end
        if is_zero_col
            push!(zero_cols, j)
        else
            push!(nonzero_cols, j)
        end
    end
    return RealisableMatroid(A, r, zero_cols, nonzero_cols, Dict{Vector{Int}, Int}())
end

function Base.convert(::Type{Matroid}, M::RealisableMatroid)
    return matroid(M)
end

@doc raw"""
    ground_set(M::RealisableMatroid)

Return the ground set of the realisable matroid `M`.
"""
function ground_set(M::RealisableMatroid)
    return Set{Int}(1:ncols(matrix(M)))
end

@doc raw"""
    flats(M::RealisableMatroid)

Return a list of all the flats of the realisable matroid `M`.
"""
function flats(M::RealisableMatroid)
    return Oscar.flats(matroid(M))
end

@doc raw"""
    rank(M::RealisableMatroid)

Return the rank of the realisable matroid `M`.
"""
function rank(M::RealisableMatroid)

    return M.rank
end

@doc raw"""
    rank(M::Matroid)

Return the rank of the matroid `M`.
"""
function rank(M::Matroid)
    return Oscar.rank(M)
end

@doc raw"""
    rank(M::RealisableMatroid, b::Set{Int})

Return the rank of the realisable matroid `M` restricted to the set `b`.
"""
function rank(M::RealisableMatroid, b::Set{Int})
    # Filter out zero columns — they cannot contribute to rank
    effective = sort!(collect(setdiff(b, M.zeroColumns)))
    isempty(effective) && return 0

    # Check cache
    cached = get(M.rankCache, effective, nothing)
    cached !== nothing && return cached

    # Compute and cache
    r = Oscar.rank(matrix(M)[:, effective])
    M.rankCache[effective] = r
    return r
end

@doc raw"""
    matroid(M::Matroid)

Return the matroid `M`.
"""
function matroid(M::Matroid)
    return M
end

@doc raw"""
    is_basis(M::RealisableMatroid, b::Set{Int})

Return `true` if `b` is a basis of the realisable matroid `M`.
"""
function is_basis(M::RealisableMatroid, b::Set{Int})
    return length(b) == rank(M)
end



###############################################################################
#
#  Analogues for Oscar matroids
#
###############################################################################

@doc raw"""
    is_basis(M::Matroid, b::Set{Int})

Return `true` if `b` is a basis of the matroid `M`.
"""
function is_basis(M::Matroid, b::Set{Int})
    return b in Set.(bases(M))
end


@doc raw"""
    ground_set(M::Union{RealisableMatroid, Matroid})

Return the ground set of a matroid as a set.
"""
function ground_set(M::Matroid)
    return Set(M.groundset)
end
