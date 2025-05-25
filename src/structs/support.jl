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

function Base.length(s::Support)
    return length(entries(s))
end

function Base.copy(s::Support)
    return Support(copy(entries(s)))
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


    # create a new dictionary to store the entries
    new_entries = Dict{Point, Height}()

    # for each point in the union, subtract the heights

    for p in points(t)
        new_entries[p] = get(entries(s), p, Nemo.PosInf()) - get(entries(t), p, Nemo.PosInf())
    end

    return Support(new_entries)

end

@doc raw"""
    update_height!(s::Support, p::Point, w::Height)

Update the height of a point in the support. Does nothing if the point is not in the support.
"""
function update_height!(s::Support, p::Point, w::Height)
    if haskey(entries(s), p)
        s.entries[p] = w
    end
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

function Base.:*(w::Height, s::Support)::Support
    return Support(Dict((p, w * h) for (p, h) in entries(s)))
end

import Base.merge

@doc raw"""
    merge(S::Support, T::Support)

Merges the support `T into the support `S`. Keeps the heights of `S` by default. Replaces points in `S` with points in `T` if they have the same entries.

This requires comparing the points by their entries.
"""
function merge(S::Support, T::Support)
    # Create a new dictionary to store the merged entries
    merged_entries = Dict{Point, Height}()

    # Copy all entries from S
    for (point, height) in S.entries
        merged_entries[point] = height
    end

    # Add entries from T, but only for points that don't exist in S
    for (pt, ht) in T.entries
        # does there exist a key with the same entries as point?
        # if there exists a key with entries(key) == entries(point), then we don't add point to the merged entries
        if !any([entries(p) == entries(pt) for p in points(S)])
            merged_entries[pt] = ht
        else
            # replace the point in T with the point in S
            for p in points(S)
                if entries(p) == entries(pt)
                    # remove pt as a key from T, and add p as a key with the height of pt
                    merged_entries[p] = S[p]
                    delete!(T.entries, pt)
                    T.entries[p] = ht
                end
            end
        end
    end

    # Return a new Support with the merged entries
    return Support(merged_entries)
end

function support(f::TropicalPolynomial)::Support

    return support(point.(collect(exponents(f))), collect(coefficients(f)))

end

import Oscar.dot
function dot(w::TropicalPoint, p::Point)
    return sum(w .* entries(p))
end

@doc raw"""
    minimum_monomials(S::Support, w::TropicalPoint)

Returns the monomials `m` in the support `S` where the minimum of $ m \cdot w + S[m] $ is achieved.
"""
function minimum_monomials(S::Support, w::TropicalPoint)
    return [m for m in points(S) if all(dot(w, m) + S[m] <= dot(w, p) + S[p] for p in points(S))]
end
