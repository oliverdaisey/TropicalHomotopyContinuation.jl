using Oscar

export Weight, convert_to_pm_type

struct Weight
    value::TropicalSemiringElem
end

function weight(w)::Weight
    return Weight(tropical_semiring()(w))
end

function Base.show(io::IO, w::Weight)
    print(io, w.value)
end

function Base.:+(a::Weight, b::Weight)::Weight
    return Weight(a.value * b.value)
end

function Base.convert(::Type{QQFieldElem}, w::Weight)
    return QQ(w.value)
end