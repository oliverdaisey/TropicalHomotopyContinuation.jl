using Oscar

export Weight, convert_to_pm_type

const Weight = Union{QQFieldElem, PosInf}

function weight(w)::Weight
    if isinf(w)
        return w 
    end
    return QQ(w)
end