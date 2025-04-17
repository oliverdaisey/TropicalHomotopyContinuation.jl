@doc raw"""
    Height

A type representing the height of a point in a tropical homotopy.
"""
const Height = Union{QQFieldElem, PosInf}

function Base.convert(::Type{Height}, w)::Height
    if isinf(w)
        return w
    end
    return QQ(w)
end
