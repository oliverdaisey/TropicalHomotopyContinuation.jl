
"""
    Tracker

A tracker for a mixed cell. Fundamental unit of the homotopy algorithm.
"""
mutable struct Tracker

    ambientSupport::MixedSupport
    activeSupport::MixedSupport
    target::MixedSupport

end

"""
    tracker(ambientSupport::MixedSupport, activeSupport::MixedSupport, target::MixedSupport)::Tracker

Construct a tracker for a mixed cell.
"""
function tracker(ambientSupport::MixedSupport, activeSupport::MixedSupport, target::MixedSupport)::Tracker
    return Tracker(ambientSupport, activeSupport, target)
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
    target(T::Tracker)

Return the target support of the tracker `T`. This is where the tracker should terminate.
"""
function target(T::Tracker)
    return T.target
end

function perturb(T::Tracker)
    # TODO: Write implementation
end

"""
    direction(T::Tracker)

Return the direction of the tracker `T`. This is the difference between the target and ambient supports.
"""
function direction(T::Tracker)
    return target(T) - ambient_support(T)
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
    weights = QQFieldElem[]
    dir = QQFieldElem[]
    for S in supports(σ)
        p1 = first(points(S))
        for p in points(S)
            if !isequal(p1, p)
                push!(rows, p1 - p)
                push!(weights, Δ[p] - Δ[p1])

                println("----------")
                println("Pushing a row corresponing to p1 = $p1 and p = $p")
                println("The row is $(p1 - p)")
                println("The difference of weights is $(Δ[p] - Δ[p1])")
                println("The direction at p is $(target(T)[p] - Δ[p])")
                println("The current weight is $(Δ[p])")
                push!(dir, target(T)[p] - Δ[p])
            end
        end
    end

    flag, inverse = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))

    if !flag
        return Nothing
    end

    return inverse * weights, inverse * dir

end