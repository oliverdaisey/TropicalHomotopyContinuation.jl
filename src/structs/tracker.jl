
"""
    Tracker

A tracker for a mixed cell. Fundamental unit of the homotopy algorithm.
"""
mutable struct Tracker

    ambientSupport::MixedSupport
    activeSupport::MixedSupport
    targets::Vector{MixedSupport}

end

"""
    tracker(ambientSupport::MixedSupport, activeSupport::MixedSupport, targets::MixedSupport)::Tracker

Construct a tracker for a mixed cell.
"""
function tracker(ambientSupport::MixedSupport, activeSupport::MixedSupport, targets::Vector{MixedSupport})::Tracker
    return Tracker(ambientSupport, activeSupport, targets)
end

"""
    active_support(T::Tracker)

Return the active support of the tracker `T`. These define the mixed cell candidate we are tracking.
"""
function active_support(T::Tracker)
    return T.activeSupport
end

"""
    ambient_support(T::Tracker)

Return the ambient support of the tracker `T`.
"""
function ambient_support(T::Tracker)
    return T.ambientSupport
end

"""
    targets(T::Tracker)

Return the targets support of the tracker `T`. This is where the tracker should terminate.
"""
function targets(T::Tracker)
    return T.targets
end

function perturb(T::Tracker)
    # TODO: Write implementation
end

"""
    direction(T::Tracker)

Return the direction of the tracker `T`. This is the difference between the targets and ambient supports.
"""
function direction(T::Tracker)
    return first(targets(T)) - ambient_support(T)
end

function Base.show(io::IO, T::Tracker)
    print(io, "Mixed cell tracker")
end

"""
    tropical_intersection_point_and_drift(T::Tracker)::Union{Nothing, Tuple{Vector{QQFieldElem}, Vector{QQFieldElem}}}

Compute the tropical intersection point and tropical drift of the tracker `T`. Returns `Nothing` if the intersection point is not well-defined.
"""
function tropical_intersection_point_and_drift(T::Tracker)::Union{Nothing, Tuple{Vector{QQFieldElem}, Vector{QQFieldElem}}}

    σ = active_support(T)
    Δ = ambient_support(T)

    rows = Vector{Int}[]
    heights = QQFieldElem[]
    dir = QQFieldElem[]
    for S in supports(σ)
        p1 = first(points(S))
        for p in points(S)
            if !isequal(p1, p)
                push!(rows, p1 - p)
                push!(heights, Δ[p] - Δ[p1])
                push!(dir, direction(T)[p] - direction(T)[p1])
            end
        end
    end

    flag, inverse = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))

    if !flag
        return Nothing
    end

    return inverse * heights, inverse * dir

end