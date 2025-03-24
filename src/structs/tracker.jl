@doc raw"""
    Tracker

A tracker for a mixed cell. Fundamental unit of the homotopy algorithm.
"""
mutable struct Tracker

    ambientSupport::MixedSupport
    mixedCells::Vector{MixedCell}
    targets::Vector{MixedSupport}

end

@doc raw"""
    tracker(ambientSupport::MixedSupport, mixedCells::Vector{MixedCellCandidate}, targets::MixedSupport)::Tracker

Construct a tracker for a mixed cell.
"""
function tracker(ambientSupport::MixedSupport, mixedCells::Vector{MixedCell}, targets::Vector{MixedSupport})::Tracker
    return Tracker(ambientSupport, mixedCells, targets)
end

@doc raw"""
    mixed_cells(T::Tracker)

Return the mixed cells of the tracker `T`.
"""
function mixed_cells(T::Tracker)
    return T.mixedCells
end

@doc raw"""
    ambient_support(T::Tracker)

Return the ambient support of the tracker `T`.
"""
function ambient_support(T::Tracker)
    return T.ambientSupport
end

@doc raw"""
    targets(T::Tracker)

Return the targets of the tracker `T`. This is a vector of mixed supports that the tracker is trying to reach.
"""
function targets(T::Tracker)
    return T.targets
end

function perturb(T::Tracker)
    # TODO: Write implementation
end

@doc raw"""
    direction(T::Tracker)

Return the direction of the tracker `T`. This is the difference between the first target and ambient support.
"""
function direction(T::Tracker)
    # TODO: If points go to infinity, make sure the points that stay finite are the same in the target and ambient support
    # TODO: In the infinity case, normalize direction to be a 0/1 vector
    return first(targets(T)) - ambient_support(T)
end

function Base.show(io::IO, T::Tracker)
    print(io, "Mixed cell tracker")
end

@doc raw"""
    tropical_intersection_point_and_drift(T::Tracker)::Union{Nothing, Tuple{Vector{QQFieldElem}, Vector{QQFieldElem}}}

Compute the tropical intersection point and tropical drift of the mixed cell σ with tracker `T`. Returns `Nothing` if the intersection point is not well-defined.
"""
function tropical_intersection_point_and_drift(T::Tracker, σ::MixedCell)::Union{Nothing,Tuple{Vector{QQFieldElem},Vector{QQFieldElem}}}

    Δ = ambient_support(T)
    τ = direction(T)
    C = chain_of_flats(σ)

    rows = Vector{Int}[]
    heights = QQFieldElem[]
    dir = QQFieldElem[]

    # add in rows from chain of flats
    pts = loopless_face(C)
    p1 = first(pts)

    for p in pts
        if !isequal(p1, p)
            push!(rows, p1 - p)
            push!(heights, 0)
            push!(dir, 0)
        end
    end

    for S in supports(σ)
        p1 = first(points(S))
        for p in points(S)
            if !isequal(p1, p)
                push!(rows, p1 - p)
                push!(heights, Δ[p] - Δ[p1])
                push!(dir, τ[p] - τ[p1])
            end
        end
    end

    flag, inverse = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))

    if !flag
        return nothing
    end

    return inverse * heights, inverse * dir
end

function tropical_intersection_point(Δ::MixedSupport)

    rows = Vector{Int}[]
    heights = QQFieldElem[]

    for S in supports(Δ)
        p1 = first(points(S))
        for p in points(S)
            if !isequal(p1, p)
                push!(rows, p1 - p)
                push!(heights, Δ[p] - Δ[p1])
            end
        end
    end

    flag, inverse = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))

    if !flag
        return nothing
    end

    return inverse * heights
end

@doc raw"""
    merge_mixed_cell!(T::Tracker, σ::MixedCell)

Merge the mixed cell `σ` into the tracker `T`. If the mixed cell is already in the tracker, do nothing.
"""
function merge_mixed_cell!(T::Tracker, σ::MixedCell)
    # check if σ is already in the tracker
    for τ in mixed_cells(T)
        # the active supports need to be the same
        if has_same_active_support(τ, σ)
            # the chains of flats need to be the same
            if chain_of_flats(τ) == chain_of_flats(σ)
                return
            end
        end
    end
    push!(T.mixedCells, σ)
end

@doc raw"""
    remove_mixed_cell!(T::Tracker, σ::MixedCell)

Remove the mixed cell `σ` from the tracker `T`.
"""
function remove_mixed_cell!(T::Tracker, σ::MixedCell)
    deleteat!(T.mixedCells, findfirst(isequal(σ), mixed_cells(T)))
end

@doc raw"""
    add_heights!(T::Tracker, Δ::MixedSupport)

Add the heights in `Δ` to the heights of the ambient support of the tracker `T`.
"""
function add_heights!(T::Tracker, Δ::MixedSupport)
    for S in supports(ambient_support(T))
        for p in points(S)
            update_height!(S, p, ambient_support(T)[p] + Δ[p])
        end
    end
end