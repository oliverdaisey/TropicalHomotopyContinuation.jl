@doc raw"""
    RealisableMatroid

An interface to a realisable matroid that provides access to the realisation matrix.
"""
struct RealisableMatroid
    realisationMatrix::MatElem{<:FieldElem}
end

function matrix(M::RealisableMatroid)
    return M.realisationMatrix
end

function matroid(M::RealisableMatroid)
    return matroid_from_matrix_columns(matrix(M))
end

function Base.show(io::IO, M::RealisableMatroid)
    print(io, "Realisable ", matroid(M))
end

function matroid(A::MatElem{<:FieldElem})
    return RealisableMatroid(A)
end

function Base.convert(::Type{Matroid}, M::RealisableMatroid)
    return matroid(M)
end

function ground_set(M::RealisableMatroid)
    return ground_set(matroid(M))
end

function flats(M::RealisableMatroid)
    return flats(matroid(M))
end

function rank(M::RealisableMatroid)
    return rank(matroid(M))
end

function matroid(M::Matroid)
    return M
end