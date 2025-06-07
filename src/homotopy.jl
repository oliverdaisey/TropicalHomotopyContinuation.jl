@doc raw"""
    move!(T::Tracker)

Move the heights of the tracker `T` to the next flip time, updating all mixed cells accordingly.

If we would reach the target before or at any flip time, do nothing.
"""
function move!(T::Tracker)

    update_number_of_moves!(T, 1)

    if isempty(mixed_cells(T))
        @vprintln :TropicalHomotopyContinuationMove 1 "No mixed cells to move. Something went wrong."
        return false
    end
    @vprintln :TropicalHomotopyContinuationMove "Moving the tracker."
    # check if we have reached the target heights
    for p in points(ambient_support(T))
        if length(targets(T)) == 1 && all(ambient_support(T)[p] == first(targets(T))[p] for p in points(first(targets(T))))
            @vprintln :TropicalHomotopyContinuationMove 1 "Reached the target heights. Tropical homotopy algorithm complete."
            return false
        end
    end

    # work out which mixed cell(s) have minimal Bergman / Jensen time
    @vprintln :TropicalHomotopyContinuationMove "Computing next flip time."
    # bergmanTimes and jensenTimes should be a dictionary mapping mixed cells to their times
    bergmanTimes = Dict{MixedCell, Height}()
    jensenTimes = Dict{MixedCell, Height}()
    for σ in mixed_cells(T)
        bergmanTimes[σ] = bergman_time(T, σ)
        jensenTimes[σ] = jensen_time(T, σ)
    end
    smallestTBergman = minimum([value for (key, value) in bergmanTimes])
    smallestTJensen = minimum([value for (key, value) in jensenTimes])
    smallestT = min(smallestTBergman, smallestTJensen) # the time at which we perform flips
    @vprintln :TropicalHomotopyContinuationMove "Next flip time is $smallestT."

    tTarget = 0

    if length(targets(T)) == 1
        tTarget = 1
    else
        tTarget = Nemo.PosInf()
    end

    if smallestT >= tTarget
        @vprintln :TropicalHomotopyContinuationMove "Flips happen after a target, so we need to check if we reached endgame or if a rebase needs to happen."
        if isinf(tTarget)
            @vprintln :TropicalHomotopyContinuationMove "Rebase needs to happen as we are sending a height to infinity."
            # set the heights to the first target heights and remove it from the list
            update_number_of_rebases!(T, 1)
            rebase!(T, first(targets(T)))
            @vprintln :TropicalHomotopyContinuationMove "Rebase complete."
            return true
        else
            @vprintln :TropicalHomotopyContinuationMove "No moves possible before reaching the target mixed subdivision. End game reached."
            # set the heights to the target heights
            add_heights!(T, direction(T))
            return false
        end
    end

    if smallestTBergman < smallestTJensen
        @vprintln :TropicalHomotopyContinuationMove "Bergman time is smaller than Jensen time."
    else
        @vprintln :TropicalHomotopyContinuationMove "Jensen time is smaller than Bergman time."
    end
    if smallestTBergman == smallestTJensen && smallestTBergman == smallestT
        update_number_of_simultaneous_bergman_and_jensen_moves!(T, 1)
        @vprintln :TropicalHomotopyContinuationMove "Both Bergman and Jensen times are equal. Performing both moves."
    end

    changingBergmanMixedCells = [σ for σ in mixed_cells(T) if bergmanTimes[σ] == smallestT]
    changingJensenMixedCells = [σ for σ in mixed_cells(T) if jensenTimes[σ] == smallestT]

    @vprintln :TropicalHomotopyContinuationMove "Mixed cells that are changing: Bergman cells: $(changingBergmanMixedCells) and Jensen cells: $(changingJensenMixedCells)."

    # these two lists should have empty intersection, if not, then a perturbation was required
    if !isempty(intersect(changingBergmanMixedCells, changingJensenMixedCells))
        @vprintln :TropicalHomotopyContinuationMove "Intersection of Bergman and Jensen mixed cells is not empty. A perturbation was required."
        perturb!(T, smallestT)
        return true
    end

    newMixedCells = MixedCell[]

    if smallestTBergman == smallestT
        @vprintln :TropicalHomotopyContinuationMove "Performing Bergman move."
        update_number_of_bergman_moves!(T, 1)
        for σ in changingBergmanMixedCells
            @vprintln :TropicalHomotopyContinuationMove "Changing mixed cell: $σ."
            bergmanFlipNewMixedCells = bergman_flip(T, σ, smallestTBergman)
            if isempty(bergmanFlipNewMixedCells)
                # this means we need to perturb
                perturb!(T, smallestT)
                return move!(T)
            end
            append!(newMixedCells, bergmanFlipNewMixedCells)
            @vprintln :TropicalHomotopyContinuationMove "$(length(newMixedCells)) new mixed cells: $newMixedCells."
            remove_mixed_cell!(T, σ)
        end
    end
    l = length(newMixedCells)
    if smallestTJensen == smallestT
        @vprintln :TropicalHomotopyContinuationMove "Performing Jensen move."
        update_number_of_jensen_moves!(T, 1)
        for σ in changingJensenMixedCells
            @vprintln :TropicalHomotopyContinuationMove "Changing mixed cell: $σ."
            jensenFlipNewMixedCells = jensen_flip(T, σ, smallestTJensen)
            if isempty(jensenFlipNewMixedCells)
                # this means we need to perturb
                perturb!(T, smallestT)
                return move!(T)
            end
            append!(newMixedCells, jensenFlipNewMixedCells)
            @vprintln :TropicalHomotopyContinuationMove "$(length(newMixedCells) - l) new mixed cells: $(newMixedCells[l+1:end])."
            remove_mixed_cell!(T, σ)
        end
    end

    # move the tracker to the heights achieved at the flip time
    @vprintln :TropicalHomotopyContinuationMove "Moving the tracker to the heights achieved at the flip time."
    add_heights!(T, smallestT*direction(T))

    # check consistency of new cells
    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuation)>0
        for σ in mixed_cells(T)
            @assert is_transverse(σ) "$(σ) is not transverse"
            @assert are_support_heights_finite(T, σ) "$(σ) has invalid mixed height data"
            @assert is_bergman_consistent(T, σ) "Cells are not Bergman consistent"
        end
    end

    @vprintln :TropicalHomotopyContinuationMove "New heights: $(dump_info(ambient_support(T)))."

    @vprintln :TropicalHomotopyContinuationMove "Merging $(length(newMixedCells)) new mixed cells."
    newMixedCellWasMerged = []
    for σ in newMixedCells
        @vprintln :TropicalHomotopyContinuationMove "Merging mixed cell: $σ."
        push!(newMixedCellWasMerged, merge_mixed_cell!(T, σ))
    end

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationMove)>1
        preFlipMultiplicities = sum(multiplicity.(changingBergmanMixedCells); init=0) + sum(multiplicity.(changingJensenMixedCells); init=0)
        postFlipMultiplicities = sum(multiplicity.(newMixedCells[findall(newMixedCellWasMerged)]); init=0)

        if preFlipMultiplicities != postFlipMultiplicities
            println("The flips in this move have changed the local multiplicities. This is not allowed.")
            println("Pre flip multiplicities: $(preFlipMultiplicities). Post flip multiplicities: $(postFlipMultiplicities).")
            println("Multiplicities of Bergman mixed cells: $(multiplicity.(changingBergmanMixedCells)).")
            println("Multiplicities of Jensen mixed cells: $(multiplicity.(changingJensenMixedCells)).")
            println("Multiplicities of new mixed cells: $(multiplicity.(newMixedCells)).")
            @assert false "Asserted, check output"
        end
    end

    # update the cached jensen and bergman times for the mixed cells that didn't change
    update_cached_times!(T, newMixedCells, smallestT)

    @vprintln :TropicalHomotopyContinuationMove "Tracker move successfully."

    update_max_mixed_cells!(T, length(mixed_cells(T)))

    return true

end

function tropical_homotopy_continuation(T::Tracker)
    while move!(T) end
end

function stable_intersection(T::Tracker) # ::Vector{TropicalPoint}
    while move!(T)
    @vprintln :TropicalHomotopyContinuation "$(T.logger)"
    @vprintln :TropicalHomotopyContinuation "$(length(mixed_cells(T))) mixed cells being tracked"
    end
    @vprintln :TropicalHomotopyContinuation "$(T.logger)"
    @vprintln :TropicalHomotopyContinuation mixed_cells(T)

    return [first(tropical_intersection_point_and_drift(T, σ)) for σ in mixed_cells(T)]

    mults = multiplicity.(mixed_cells(T))
    pts = [first(tropical_intersection_point_and_drift(T, σ)) for σ in mixed_cells(T)]

    uniquePts = unique(pts)
    uniqueMults = zeros(Int, length(uniquePts))

    for (p, m) in zip(pts, mults)
        i = findfirst(t -> is_equal(p, t), uniquePts)
        uniqueMults[i] += m
    end

    return uniquePts, uniqueMults
end
