
"""
    DualSupport{dualType}(points::Matrix{Int})

A support for a dual object, with rows as points.

"""
mutable struct DualSupport{dualType<:DualType}
    points::Matrix{Int}

    function DualSupport{dualType}(points::Matrix{Int}) where {dualType<:DualType}
        return new{dualType}(points)
    end

    function DualSupport{Hypersurface}(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}})
        return new{Hypersurface}(generate_support(f))
    end
end

mutable struct LazyDualSupport{dualType<:Union{Linear, InvertedLinear}}

    # in this case, points are not explicitly given, but there should be a way to generate them
    n::Int
    k::Int

    function LazyDualSupport{dualType}(n::Int, k::Int) where {dualType<:Union{Linear, InvertedLinear}}
        return new{dualType}(n, k)
    end
end

function points(s::DualSupport)
    return s.points
end

function tropical_ambient_dim(s::DualSupport)
    return size(s.points, 2)
end

function Base.getindex(s::DualSupport, i::Int)
    return s.points[i, :]
end

function Base.getindex(s::DualSupport, i::Int, j::Int)
    return s.points[i, j]
end

function Base.getindex(s::DualSupport, I::Vector{Int})
    return s.points[I, :]
end

function Base.getindex(s::DualSupport, i::Int, ::Colon)
    return s.points[i, :]
end

function Base.size(s::DualSupport, i::Int)
    return size(s.points, i)
end

function Base.length(s::DualSupport)
    return size(s.points, 1)
end

function Base.iterate(s::DualSupport)
    return iterate(s[i] for i in 1:length(s))
end

function Base.show(io::IO, s::DualSupport)
    print(io, "Dual support of type $(typeof(s).parameters[1]) with points $(s.points)")
end

function Base.:(==)(s1::DualSupport, s2::DualSupport)
    return s1.points == s2.points && typeof(s1) == typeof(s2)
end

"""
    generate_support(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}})

Generate the support of a tropical polynomial.

# Arguments
- `f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}`: The tropical polynomial.

# Returns
- The support of the tropical polynomial.
"""
function generate_support(f::AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}})::Matrix{Int64}
    return transpose(hcat(collect(exponents(f))...))
end

function tropical_codim(S::DualSupport{Hypersurface})
    return 1
end

function tropical_codim(S::DualSupport{<:Union{Linear, InvertedLinear}})
    return length(findall(x -> x != 0, S[1]))
end