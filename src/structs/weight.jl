using Oscar

export Weight

struct Weight
    value::TropicalSemiringElem
end

function weight(w)::Weight
    return Weight(tropical_semiring()(w))
end

function Base.show(io::IO, w::Weight)
    print(io, w.value)
end