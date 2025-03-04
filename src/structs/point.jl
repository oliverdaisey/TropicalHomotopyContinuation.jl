export Point, point, minkowski_sum

@doc raw"""
    Point

A point in the integer lattice. This is a simple wrapper around a vector of integers. Encodes exponent vectors and monomials.
"""
struct Point
    entries::AbstractVector{Int}
end

@doc raw"""
    point(args...)

Construct a point from a list of integers.
"""
function point(args...)::Point
    return Point(collect(args))
end

@doc raw"""
    point(arg::AbstractVector{Int})

Construct a point from a vector of integers.
"""
function point(arg::AbstractVector{Int})
    return Point(arg)
end

function Base.show(io::IO, p::Point)
    print(io, "($(join(p.entries, ", ")))")
end

Base.length(p::Point) =length(p.entries)

@doc raw"""
    convert(Vector, p::Point)

Convert a point to a vector.
"""
function Base.convert(::Any, p::Point)
    return collect(p.entries)
end

function Base.getindex(p::Point, i::Int)
    return p.entries[i]
end

function Base.setindex!(p::Point, v, i::Int)
    p.entries[i] = v
end

function Base.:+(a::Point, b::Point)::Point
    return Point(a.entries .+ b.entries)
end

function Base.:-(a::Point, b::Point)::Point
    return Point(a.entries .- b.entries)
end

@doc raw"""
    minkowski_sum(args...)

Compute the Minkowski sum of a list of vector-like objects.
"""
function minkowski_sum(args::Vector...)
    return vcat([sum(vecs) for vecs in Iterators.product(args...)]...)
end