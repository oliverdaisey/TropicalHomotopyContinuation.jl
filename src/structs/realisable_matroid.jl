@doc raw"""
    RealisableMatroid

An interface to a realisable matroid that provides access to the realisation matrix.
"""
struct RealisableMatroid
    realisationMatrix::MatElem{<:FieldElem}
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

function Base.show(io::IO, M::RealisableMatroid)
    print(io, "Realisable ", matroid(M))
end

@doc raw"""
    matroid(A::MatElem{<:FieldElem})

Construct a realisable matroid from the realisation matrix `A`.
"""
function matroid(A::MatElem{<:FieldElem})
    return RealisableMatroid(A)
end

function Base.convert(::Type{Matroid}, M::RealisableMatroid)
    return matroid(M)
end

@doc raw"""
    ground_set(M::RealisableMatroid)

Return the ground set of the realisable matroid `M`.
"""
function ground_set(M::RealisableMatroid)
    return ground_set(matroid(M))
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
    return Oscar.rank(matroid(M))
end

@doc raw"""
    rank(M::Matroid)

Return the rank of the matroid `M`.
"""
function rank(M::Matroid)
    return Oscar.rank(M)
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
    return isequal(Oscar.rank(matrix(M)[:, collect(b)]), Oscar.rank(matrix(M)))
end

@doc raw"""
    is_basis(M::Matroid, b::Set{Int})

Return `true` if `b` is a basis of the matroid `M`.
"""
function is_basis(M::Matroid, b::Set{Int})
    return b in Set.(bases(M))
end