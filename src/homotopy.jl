@doc raw"""
    move!(T::Tracker)

Move the heights of the tracker `T` to the next flip time, updating all mixed cells accordingly.

If we would reach the target before or at any flip time, do nothing.
"""
function move!(T::Tracker)
    
    @debug "Moving the tracker."
    # check if we have reached the target heights
    for p in points(ambient_support(T))
        if length(targets(T)) == 1 && all(ambient_support(T)[p] == first(targets(T))[p] for p in points(first(targets(T))))
            @debug "Reached the target heights. Tropical homotopy algorithm complete."
            return
        end
    end

    # work out which mixed cell(s) have minimal Bergman / Jensen time
    @debug "Computing next flip time."
    smallestTBergman = minimum([bergman_time(T, σ) for σ in mixed_cells(T)])
    smallestTJensen = minimum([jensen_time(T, σ) for σ in mixed_cells(T)])
    smallestT = min(smallestTBergman, smallestTJensen) # the time at which we perform flips
    @debug "Next flip time is $smallestT."

    tTarget = 0

    if length(targets(T)) == 1
        tTarget = 1
    else
        tTarget = Nemo.PosInf()
    end

    if smallestT >= tTarget
        @debug "Flips happen after a target, so we need to check if we reached endgame or if a rebase needs to happen."
        if length(targets(T)) != 1
            @debug "Rebase needs to happen as we are sending a height to infinity."
            # set the heights to the first target heights and remove it from the list
            rebase!(T, first(targets(T)))
            @debug "Rebase complete."
            return true
        else
            @debug "No moves possible before reaching the target mixed subdivision. End game reached."
            # set the heights to the target heights
            add_heights!(T, direction(T))
            return false
        end
    end

    if smallestTBergman < smallestTJensen
        @debug "Bergman time is smaller than Jensen time."
    else
        @debug "Jensen time is smaller than Bergman time."
    end
    if smallestTBergman == smallestTJensen && smallestTBergman == smallestT 
        @debug "Both Bergman and Jensen times are equal. Performing both moves."
    end

    changingBergmanMixedCells = [σ for σ in mixed_cells(T) if bergman_time(T, σ) == smallestT]
    changingJensenMixedCells = [σ for σ in mixed_cells(T) if jensen_time(T, σ) == smallestT]

    @debug "Mixed cells that are changing: Bergman cells: $(changingBergmanMixedCells) and Jensen cells: $(changingJensenMixedCells)."

    # these two lists should have empty intersection, if not, then a perturbation was required
    @assert isempty(intersect(changingBergmanMixedCells, changingJensenMixedCells)) "A perturbation is required, but these are not implemented yet. Violating mixed cell(s): $(intersect(changingBergmanMixedCells, changingJensenMixedCells)). Bergman time(s): $(bergman_time.(Ref(T), intersect(changingBergmanMixedCells, changingJensenMixedCells))). Jensen time(s): $(jensen_time.(Ref(T), intersect(changingBergmanMixedCells, changingJensenMixedCells)))"

    newMixedCells = MixedCell[]

    if smallestTBergman <= smallestTJensen
        @debug "Performing Bergman moves."
        for σ in changingBergmanMixedCells
            @debug "Changing mixed cell: $σ."
            push!(newMixedCells, bergman_flip(T, σ)...)
            @debug "$(length(newMixedCells)) new mixed cells: $newMixedCells."
            remove_mixed_cell!(T, σ)
        end
    end
    l = length(newMixedCells)
    if smallestTJensen <= smallestTBergman
        @debug "Performing Jensen moves."
        for σ in changingJensenMixedCells
            @debug "Changing mixed cell: $σ."
            push!(newMixedCells, jensen_flip(T, σ)...)
            @debug "$(length(newMixedCells) - l) new mixed cells: $(newMixedCells[l+1:end])."
            remove_mixed_cell!(T, σ)
        end
    end

    # move the tracker to the heights achieved at the flip time
    @debug "Moving the tracker to the heights achieved at the flip time."
    add_heights!(T, min(smallestTBergman, smallestTJensen)*direction(T))

    @debug "New heights: $(dump_info(ambient_support(T)))."

    @debug "Merging $(length(newMixedCells)) new mixed cells."
    for σ in newMixedCells
        @debug "Merging mixed cell: $σ."
        merge_mixed_cell!(T, σ)
    end

    @debug "Tracker move successfully."

    return true

end

function tropical_homotopy_continuation(T::Tracker)
    while move!(T) end
end

function stable_intersection(T::Tracker)::Vector{TropicalPoint}
    while move!(T) end
    return [first(tropical_intersection_point_and_drift(T, σ)) for σ in mixed_cells(T)]
end