export vector_of_points

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
    return keys(entries(s))
end

function weights(s::Support)
    return values(entries(s))
end

function Base.show(io::IO, s::Support)

    print(io, "Dual support $(join(["$(p)↑$(w)" for (p, w) in s.entries], ", "))")
end

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