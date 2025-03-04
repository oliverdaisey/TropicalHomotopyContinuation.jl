using Oscar

export Height, height

const Height = Union{QQFieldElem, PosInf}

function height(w)::Height
    if isinf(w)
        return w 
    end
    return QQ(w)
end