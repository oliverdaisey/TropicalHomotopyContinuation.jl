using Oscar

export Weight, convert_to_pm_type

"""
    Weight

A weight is an element of the tropical semiring. This is used to define the height of a point in a support.
"""
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
    return Weight(a.value + b.value)
end

function Base.:-(a::Weight, b::Weight)::Weight
    return Weight(a.value / b.value)
end

function Base.:*(a::Weight, b::Weight)::Weight
    return Weight(a.value * b.value)
end

function Base.convert(::Type{QQFieldElem}, w::Weight)
    return QQ(w.value)
end

function Base.isless(a::Weight, b)
    return a.value < tropical_semiring()(b)
end

function Base.isequal(a::Weight, b)
    return a.value == tropical_semiring()(b)
end

function convert_to_pm_type(w::Weight)
    return tropical_semiring()(w.value)
end