mutable struct Support

    entries::Dict{Point, Weight}

end

function support(points, weights)::Support
    @assert length(points) == length(weights) "The number of points and weights must be the same"
    @assert all(p -> length(p) == length(points[1]), points) "All points must have the same dimension"

    return Support(Dict(zip(points, weights)))
end

function points(s::Support)
    return keys(s.entries)
end

function weights(s::Support)
    return values(s.entries)
end

function Base.show(io::IO, s::Support)
    print(io, "Dual support with points $(join(points(s), ", ")) and weights $(join(weights(s), ", "))")
end

function update_weight(s::Support, p::Point, w::Weight)
    @assert haskey(s.entries, p) "The point $p is not in the support"
    s.entries[p] = w
end