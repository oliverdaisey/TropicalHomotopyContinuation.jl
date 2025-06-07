@doc raw"""
    Tracker

A tracker for a mixed cell. Fundamental unit of the homotopy algorithm.
"""
mutable struct Tracker

    ambientSupport::MixedSupport
    mixedCells::Vector{MixedCell}
    bergmanTimes::Dict{MixedCell, Height}
    jensenTimes::Dict{MixedCell, Height}
    targets::Vector{MixedSupport}
    logger::Logger # logs data about the tracker

end

@doc raw"""
    tracker(ambientSupport::MixedSupport, mixedCells::Vector{MixedCellCandidate}, targets::MixedSupport)::Tracker

Construct a tracker for a mixed cell.
"""
function tracker(ambientSupport::MixedSupport, mixedCells::Vector{MixedCell}, targets::Vector{MixedSupport})::Tracker
    T = Tracker(copy(ambientSupport), copy(mixedCells), Dict{MixedCell, Height}(), Dict{MixedCell, Height}(), copy(targets), logger())
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

function perturb!(T::Tracker, timeOfFailure::Height)

    @assert !isinf(timeOfFailure) "Time of failure should be finite"
    @vprintln :TropicalHomotopyContinuationPerturb "time of failure = $(timeOfFailure)"
    dir = direction(T)
    add_heights!(T, (timeOfFailure / 2)*dir)
    # choose a random point in height space
    targetSupports = tuple(copy.(supports(ambient_support(T)))...)
    # update all the heights to random values
    for S in targetSupports
        for p in sort(points(S))
            if !isinf(S[p])
                newHeight = QQ(rand(Int16))
                update_height!(S, p, newHeight)
            end
        end
    end
    targetMixedSupport = mixed_support(targetSupports)


    newTracker = tracker(ambient_support(T), mixed_cells(T), [targetMixedSupport])

    bergmanTimes = Dict{MixedCell, Height}()
    jensenTimes = Dict{MixedCell, Height}()

    check_cached_times(newTracker)

    for σ in mixed_cells(newTracker)
        bergmanTimes[σ] = compute_bergman_time(newTracker, σ)
        jensenTimes[σ] = compute_jensen_time(newTracker, σ)
    end

    check_cached_times(newTracker)

    smallestTBergman = minimum([value for (key, value) in bergmanTimes])
    smallestTJensen = minimum([value for (key, value) in jensenTimes])
    smallestT = min(smallestTBergman, smallestTJensen) # the time at which we perform flips
    @vprintln :TropicalHomotopyContinuationPerturb "smallestT = $(smallestT)"
    @vprintln :TropicalHomotopyContinuationPerturb "smallestTBergman = $(smallestTBergman)"
    @vprintln :TropicalHomotopyContinuationPerturb "smallestTJensen = $(smallestTJensen)"

    dir = direction(newTracker)
    @v_do :TropicalHomotopyContinuationPerturb show_heights(T)
    @vprintln :TropicalHomotopyContinuationPerturb "Showing heights in direction between heights after first rebase and targetMixedSupport"
    @v_do :TropicalHomotopyContinuationPerturb show_heights(dir)
    add_heights!(T, (smallestT/2)*dir)
    @v_do :TropicalHomotopyContinuationPerturb show_heights(T)
    add_heights!(first(targets(T)), (smallestT/2)*dir)

    # delete all cached times
    for σ in mixed_cells(T)
        delete!(T.bergmanTimes, σ)
        delete!(T.jensenTimes, σ)
    end

end

@doc raw"""
    direction(T::Tracker)

Return the direction of the tracker `T`. This is the difference between the first target and ambient support.
"""
function direction(T::Tracker)
    dir = first(targets(T)) - ambient_support(T)
    if any(isinf.(heights(dir)))
        # every height should be either infinity or 0
        @assert all(h -> h == 0 || isinf(h), heights(dir)) "Invalid direction for tracker (the direction is $(supports(dir)))"
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

function Base.:+(::Nemo.PosInf, ::Nemo.PosInf)
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
        if !is_equal(p1, p)
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
            if !is_equal(p1, p) # && !isinf(Δ[p])
                push!(rows, p1 - p)
                push!(heights, Δ[p] - Δ[p1])
                push!(dir, τ[p] - τ[p1])
            end
        end
    end

    flag, inverse = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))

    @assert flag "Matrix not invertible, $(σ) positive-dimensional"

    return inverse * heights, inverse * dir
end


@doc raw"""
    merge_mixed_cell!(T::Tracker, σ::MixedCell)

Merge the mixed cell `σ` into the tracker `T`. If the mixed cell is already in the tracker, do nothing.

Returns `true` if the mixed cell was added, `false` otherwise.
"""
function merge_mixed_cell!(T::Tracker, σ::MixedCell)::Bool
    # check if σ is already in the tracker
    for τ in mixed_cells(T)
        # the active supports need to be the same
        if has_same_active_support(τ, σ)
            # the chains of flats need to be the same
            if is_equal(chain_of_flats(τ), chain_of_flats(σ))
                return false
            end
        end
    end
    push!(T.mixedCells, σ)
    return true
end

@doc raw"""
    remove_mixed_cell!(T::Tracker, σ::MixedCell)

Remove the mixed cell `σ` from the tracker `T`.
"""
function remove_mixed_cell!(T::Tracker, σ::MixedCell)
    T.mixedCells = setdiff(T.mixedCells, [σ])
    # remove the mixed cell from the bergman and jensen times
    delete!(T.bergmanTimes, σ)
    delete!(T.jensenTimes, σ)
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

function add_heights!(Δ1::MixedSupport, Δ2::MixedSupport)
    for S in supports(Δ1)
        for p in points(S)
            update_height!(S, p, Δ1[p] + Δ2[p])
        end
    end
end

@doc raw"""
    tracker(startingSupport::MixedSupport, targetSupport::MixedSupport, mixedCells::Vector{MixedCell})::Tracker

Construct a tracker for a list of mixed cells. Tracks heights that need to go to infinity one at a time. Then tracks finite heights along a straight line to the target heights all at once.

The starting support needs to be a superset of the target support.
"""
function tracker(startingSupport::MixedSupport, targetSupport::MixedSupport, mixedCells::Vector{MixedCell}; path::Symbol=:coefficient_wise)::Tracker

    @assert is_subset(targetSupport, startingSupport) "The target support needs to be a subset of the starting support."

    if path == :coefficient_wise
        targets = coefficientwise_homotopy(startingSupport, targetSupport)
    elseif path == :straight_line
        targets = straight_line_homotopy(startingSupport, targetSupport)
    else
        error("invalid path strategy")
    end

    return tracker(startingSupport, mixedCells, targets)

end

function coefficientwise_homotopy(startingSupport::MixedSupport, targetSupport::MixedSupport)::Vector{MixedSupport}
    targets = Vector{MixedSupport}()

    # temporarily add the starting support to work out which heights need to go to infinity
    push!(targets, startingSupport)
    for p in points(startingSupport)
        if isinf(targetSupport[p])
            # take the latest support in `targets` and set the height of `p` to infinity
            newTarget = copy(last(targets))
            update_height!(newTarget, p, Nemo.PosInf())
            push!(targets, newTarget)
        end
    end

    # at this point no more heights need to be set to infinity, so go to the target heights all at once
    # TODO: this should probably also be coefficientwise
    push!(targets, targetSupport)

    # do not return the first support in the list, it is the starting support
    return targets[2:end]
end

function straight_line_homotopy(startingSupport::MixedSupport, targetSupport::MixedSupport)::Vector{MixedSupport}
    targets = Vector{MixedSupport}()

    # stage 1: send heights of points in startingSupport but not in targetSupport to infinity
    targetMidway = copy(startingSupport)
    for p in points(startingSupport)
        if isinf(targetSupport[p])
            update_height!(targetMidway, p, Nemo.PosInf())
        end
    end
    push!(targets, targetMidway)

    # stage 2: send heights of points in startingSupport and targetSupport to the correct value
    push!(targets, targetSupport)

    return targets

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

    # invalidate the cached times
    for σ in mixed_cells(T)
        delete!(T.bergmanTimes, σ)
        delete!(T.jensenTimes, σ)
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

function update_bergman_time!(T::Tracker, σ::MixedCell, bergmanTime::Height)
    T.bergmanTimes[σ] = bergmanTime
end

function update_jensen_time!(T::Tracker, σ::MixedCell, jensenTime::Height)
    T.jensenTimes[σ] = jensenTime
end

function bergman_time(T::Tracker, σ::MixedCell)::Height
    if haskey(T.bergmanTimes, σ)
        return T.bergmanTimes[σ]
    else
        return compute_bergman_time(T, σ)
    end
end

function jensen_time(T::Tracker, σ::MixedCell)::Height
    if haskey(T.jensenTimes, σ)
        return T.jensenTimes[σ]
    else
        return compute_jensen_time(T, σ)
    end
end

@doc raw"""
    update_cached_times!(T::Tracker, newMixedCells::Vector{MixedCell}, smallestT::Height)

Update the cached jensen and bergman times for the mixed cells that didn't change.
"""
function update_cached_times!(T::Tracker, newMixedCells::Vector{MixedCell}, smallestT::Height)
    # update the cached jensen and bergman times for the mixed cells that didn't change
    jensenTimes = T.jensenTimes
    bergmanTimes = T.bergmanTimes

    for σ in mixed_cells(T)
        if σ in newMixedCells
            continue
        end
        # if we get to this point, then σ is not in newMixedCells
        # if we already cached a time then update it
        if length(targets(T)) == 1
            # in this case times are defined as the fraction of the distance between the current heights and target heights
            # so updating the cache should be done multiplicatively and not additively
            if haskey(jensenTimes, σ) && !isinf(jensenTimes[σ])
                update_jensen_time!(T, σ, (jensenTimes[σ] - smallestT) / (1 - smallestT))
            end
            if haskey(bergmanTimes, σ) && !isinf(bergmanTimes[σ])
                update_bergman_time!(T, σ, (bergmanTimes[σ] - smallestT) / (1 - smallestT))
            end
            continue
        end
        if haskey(jensenTimes, σ) && !isinf(jensenTimes[σ])
            update_jensen_time!(T, σ, jensenTimes[σ] - smallestT)
        end
        if haskey(bergmanTimes, σ) && !isinf(bergmanTimes[σ])
            update_bergman_time!(T, σ, bergmanTimes[σ] - smallestT)
        end
    end

end

function show_heights(T::Tracker)
    println("----------------------")
    for p in sort(points(ambient_support(T)))
        println("$(p) has height $(ambient_support(T)[p])")
    end
    println("----------------------")
end

function show_heights(Δ::MixedSupport)
    println("----------------------")
    for p in sort(points(Δ))
        println("$(p) has height $(Δ[p])")
    end
    println("----------------------")
end

function check_cached_times(T::Tracker)
    for σ in mixed_cells(T)
        if haskey(T.bergmanTimes, σ)
            @assert T.bergmanTimes[σ] == compute_bergman_time(T, σ, disable_cache=true) "Bergman time for $(σ) is not cached correctly"
        end
        if haskey(T.jensenTimes, σ)
            @assert T.jensenTimes[σ] == compute_jensen_time(T, σ, disable_cache=true) "Jensen time for $(σ) is not cached correctly"
        end
    end
end
