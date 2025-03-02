export vector_of_points, points

"""
    Support

A support is a collection of points with associated weights.
"""
struct Support

    entries::Dict{Point, Weight}

end

function entries(s::Support)
    return s.entries
end

function support(points, weights)::Support
    @assert length(points) == length(weights) "The number of points and weights must be the same"
    @assert all(p -> length(p) == length(points[1]), points) "All points must have the same dimension"

    return Support(Dict(zip(points, weights)))
end

function points(s::Support)
    return collect(keys(entries(s)))
end

function weights(s::Support)
    return collect(values(entries(s)))
end

function Base.show(io::IO, s::Support)

    print(io, "Dual support $(join(["$(p)↑$(w)" for (p, w) in s.entries], ", "))")
end

"""
    update_weight!(s::Support, p::Point, w::Weight)

Update the weight of a point in the support.
"""
function update_weight!(s::Support, p::Point, w::Weight)
    @assert haskey(entries(s), p) "The point $p is not in the support"
    s.entries[p] = w
end

"""
    vector_of_points(s::Support)

Convenience function to return the points of a support as a vector.
"""
function vector_of_points(s::Support)
    return collect(points(s))
end

function Base.in(p::Point, s::Support)
    return haskey(entries(s), p)
end

"""
    is_subset(S::Support, T::Support)

Returns true if the support S is a subset of the support T.
"""
function is_subset(S::Support, T::Support)
   return all(p -> in(p, T), points(S))
end

"""
    support(S::Support, p::Point)

Merges the point `p` into the support `S`.
"""
function support(S::Support, p::Point)::Support
    return Support(merge(Dict(p => weight(1)), entries(S)))
end

function Base.getindex(s::Support, p::Point)
    return get(entries(s), p, tropical_semiring()(one(QQ)))
end