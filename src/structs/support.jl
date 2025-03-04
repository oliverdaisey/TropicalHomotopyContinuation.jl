export vector_of_points, points

@doc raw"""
    Support

A support is a collection of points with associated heights.
"""
struct Support

    entries::Dict{Point, Height}

end

function entries(s::Support)
    return s.entries
end

function support(points, heights)::Support
    @assert length(points) == length(heights) "The number of points and heights must be the same"
    @assert all(p -> length(p) == length(points[1]), points) "All points must have the same dimension"

    return Support(Dict(zip(points, heights)))
end

function points(s::Support)
    return collect(keys(entries(s)))
end

function heights(s::Support)
    return collect(values(entries(s)))
end

function Base.show(io::IO, s::Support)

    print(io, "Dual support $(join(["$(p)↑$(w)" for (p, w) in s.entries], ", "))")
end

function Base.:-(s::Support, t::Support)::Support

    @assert length(points(s)) == length(points(t)) "The number of points in the supports must be the same"
    return Support(mergewith(-)(entries(s), entries(t)))
    
end

@doc raw"""
    update_height!(s::Support, p::Point, w::Height)

Update the height of a point in the support.
"""
function update_height!(s::Support, p::Point, w::Height)
    @assert haskey(entries(s), p) "The point $p is not in the support"
    s.entries[p] = w
end

@doc raw"""
    vector_of_points(s::Support)

Convenience function to return the points of a support as a vector.
"""
function vector_of_points(s::Support)
    return collect(points(s))
end

function Base.in(p::Point, s::Support)
    return haskey(entries(s), p)
end

@doc raw"""
    is_subset(S::Support, T::Support)

Returns true if the support S is a subset of the support T.
"""
function is_subset(S::Support, T::Support)
   return all(p -> in(p, T), points(S))
end

@doc raw"""
    support(S::Support, p::Point)

Merges the point `p` into the support `S`. Gives height 0 by default.
"""
function support(S::Support, p::Point)::Support
    return Support(merge(Dict(p => 0), entries(S)))
end

function Base.getindex(s::Support, p::Point)
    return get(entries(s), p, tropical_semiring()(one(QQ)))
end