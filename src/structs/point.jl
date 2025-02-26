export Point

struct Point
    entries::AbstractVector{Int}
end

function point(args...)::Point
    return Point(collect(args))
end

function Base.show(io::IO, p::Point)
    print(io, "($(join(p.entries, ", ")))")
end

function Base.length(p::Point)
    return length(p.entries)
end