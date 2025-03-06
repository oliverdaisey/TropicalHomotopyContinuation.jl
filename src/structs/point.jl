export Point, point, minkowski_sum, entries

const TropicalPoint = Union{Vector{QQFieldElem}, Vector{Int}}

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

function Base.convert(::Type{Vector{QQFieldElem}}, v::Vector{Int})
    return QQ.(v)
end

function Base.convert(::Type{Point}, v::Vector{QQFieldElem})
    return point(convert(Vector{Int}, v))
end

function entries(p::Point)
    return p.entries
end

@doc raw"""
    matrix(v::Vector{Point})

Convert a vector with elements of type `Point` to an Oscar matrix.
"""
function matrix(v::Vector{Point})
    return Oscar.matrix(QQ, [entries(p) for p in v])
end