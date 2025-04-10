export Tracker, tracker, mixed_cells, ambient_support, targets, direction, tropical_intersection_point_and_drift, tropical_intersection_point

@doc raw"""
    Tracker

A tracker for a mixed cell. Fundamental unit of the homotopy algorithm.
"""
mutable struct Tracker

    ambientSupport::MixedSupport
    mixedCells::Vector{MixedCell}
    targets::Vector{MixedSupport}
    logger::Logger # logs data about the tracker

end

@doc raw"""
    tracker(ambientSupport::MixedSupport, mixedCells::Vector{MixedCellCandidate}, targets::MixedSupport)::Tracker

Construct a tracker for a mixed cell.
"""
function tracker(ambientSupport::MixedSupport, mixedCells::Vector{MixedCell}, targets::Vector{MixedSupport})::Tracker
    T = Tracker(ambientSupport, mixedCells, targets, logger())
    update_max_mixed_cells!(T, length(mixedCells))
    return T
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
    dir = first(targets(T)) - ambient_support(T)
    if any(isinf.(heights(dir)))
        # every height should be either infinity or 0
        @assert all(h -> h == 0 || isinf(h), heights(dir)) "Invalid direction for tracker"
        for p in points(dir)
            if isinf(dir[p])
                update_height!(dir, p, QQ(1))
            end
        end
    end

    return dir
end

function Base.:-(::Nemo.PosInf, ::Nemo.QQFieldElem)
    return Nemo.PosInf()
end

function Base.show(io::IO, T::Tracker)
    print(io, "Mixed cell tracker")
end

@doc raw"""
    tropical_intersection_point_and_drift(T::Tracker)::Union{Nothing, Tuple{Vector{QQFieldElem}, Vector{QQFieldElem}}}

Compute the tropical intersection point and tropical drift of the mixed cell σ with tracker `T`.
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
        # get the first point in S that has nonzero height
        p1 = first(points(S))
        for p in points(S)
            if !isinf(Δ[p])
                p1 = p
                break
            end
        end

        for p in points(S)
            if !isequal(p1, p) # && !isinf(Δ[p])
                push!(rows, p1 - p)
                push!(heights, Δ[p] - Δ[p1])
                push!(dir, τ[p] - τ[p1])
            end
        end
    end

    flag, inverse = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))
    
    if !flag
        println(σ, " is not a valid mixed cell")
    end
    @assert flag "Matrix not invertible"

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
            if isequal(chain_of_flats(τ), chain_of_flats(σ))
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
    T.mixedCells = setdiff(T.mixedCells, [σ])
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

@doc raw"""
    tracker(startingSupport::MixedSupport, targetSupport::MixedSupport, mixedCells::Vector{MixedCell})::Tracker

Construct a tracker for a list of mixed cells. Tracks heights that need to go to infinity one at a time. Then tracks finite heights along a straight line to the target heights all at once.

The starting support needs to be a superset of the target support.
"""
function tracker(startingSupport::MixedSupport, targetSupport::MixedSupport, mixedCells::Vector{MixedCell})::Tracker

    @assert is_subset(targetSupport, startingSupport) "The target support needs to be a subset of the starting support."

    targets = Vector{MixedSupport}()
    push!(targets, startingSupport)
    # Work out which heights need to go to infinity
    for p in points(startingSupport)
        if isinf(targetSupport[p])
            # take the latest support in `targets` and set the height of `p` to infinity
            newTarget = copy(last(targets))
            update_height!(newTarget, p, Nemo.PosInf())
            push!(targets, newTarget)
        end
    end

    # at this point no more heights need to be set to infinity, so go to the target heights all at once
    push!(targets, targetSupport)

    # the first support in the list is the starting support, remove it
    targets = targets[2:end]

    return tracker(startingSupport, mixedCells, targets)

    
end

function rebase!(T::Tracker, Δ::MixedSupport)
    # remove any points from T whose height is infinity in Δ

    for p in points(ambient_support(T))
        if isinf(Δ[p])
            for S in supports(ambient_support(T))
                if in(p, S)
                    # remove p from the keys of the entries of S
                    delete!(entries(S), p)
                end
            end

            # if any of the mixed cells involve p, remove them

            for σ in mixed_cells(T)
                if p in points(σ)
                    update_number_of_diverged_mixed_cells!(T, 1)
                    remove_mixed_cell!(T, σ)
                end
            end
        end
    end

    # remove the first target from the list of targets
    T.targets = T.targets[2:end]
    # println("Number of targets left: ", length(T.targets))

    for σ in mixed_cells(T)
        # check that the matrix coming from σ is invertible
        @assert is_transverse(σ) "$(σ) is not transverse"
        @assert are_support_heights_finite(T, σ) "$(σ) has invalid mixed height data"
    end
end

function update_max_mixed_cells!(T::Tracker, n::Int)
    T.logger.maxMixedCells = max(T.logger.maxMixedCells, n)
end

function update_number_of_jensen_moves!(T::Tracker, n::Int)
    T.logger.numberOfJensenMoves += n
end
function update_number_of_bergman_moves!(T::Tracker, n::Int)
    T.logger.numberOfBergmanMoves += n
end

function update_number_of_diverged_mixed_cells!(T::Tracker, n::Int)
    T.logger.numberOfDivergedMixedCells += n
end

function update_number_of_moves!(T::Tracker, n::Int)
    T.logger.numberOfMoves += n
end

function update_number_of_rebases!(T::Tracker, n::Int)
    T.logger.numberOfRebases += n
end

function update_number_of_simultaneous_bergman_and_jensen_moves!(T::Tracker, n::Int)
    T.logger.numberOfSimultaneousBergmanAndJensenMoves += n
end

function are_support_heights_finite(T::Tracker, σ::MixedCell)
    Δ = ambient_support(T)
    for p in points(active_support(σ))
        if isinf(Δ[p])
            return false
        end
    end
    
    return true

end

function is_bergman_consistent(T::Tracker, σ::MixedCell)
    
    chainOfFlats = chain_of_flats(σ)
    w, _ = tropical_intersection_point_and_drift(T, σ)
    inequalities, equalities = cone(chainOfFlats)

    # check that we are inside the cone
    if length(equalities) > 0
        @assert all([sum(w .* v) == 0 for v in equalities]) "The intersection point is not in the cone (equality violated)"
        @assert all([sum(w .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone (inequality violated) intersection point: $(w) chain of flats: $(chainOfFlats)"
    else
        @assert all([sum(w .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone"
    end

    return true

end