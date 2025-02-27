export Point, point, minkowski_sum

struct Point
    entries::AbstractVector{Int}
end

function point(args...)::Point
    return Point(collect(args))
end

function point(arg::AbstractVector{Int})
    return Point(arg)
end

function Base.show(io::IO, p::Point)
    print(io, "($(join(p.entries, ", ")))")
end

function Base.length(p::Point)
    return length(p.entries)
end

function Base.convert(Vector, p::Point)
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

function minkowski_sum(args::Vector...)
    return vcat([sum(vecs) for vecs in Iterators.product(args...)]...)
end