using Oscar

export Height

const Height = Union{QQFieldElem, PosInf}

function Base.convert(::Type{Height}, w)::Height
    if isinf(w)
        return w 
    end
    return QQ(w)
end