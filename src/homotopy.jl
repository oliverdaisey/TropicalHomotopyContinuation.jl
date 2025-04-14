export move!, tropical_homotopy_continuation, stable_intersection

@doc raw"""
    move!(T::Tracker)

Move the heights of the tracker `T` to the next flip time, updating all mixed cells accordingly.

If we would reach the target before or at any flip time, do nothing.
"""
function move!(T::Tracker)
    
    update_number_of_moves!(T, 1)
    
    if isempty(mixed_cells(T))
        @error "No mixed cells to move. Something went wrong."
        return false
    end
    @debug "Moving the tracker."
    # check if we have reached the target heights
    for p in points(ambient_support(T))
        if length(targets(T)) == 1 && all(ambient_support(T)[p] == first(targets(T))[p] for p in points(first(targets(T))))
            @debug "Reached the target heights. Tropical homotopy algorithm complete."
            return false
        end
    end

    # work out which mixed cell(s) have minimal Bergman / Jensen time
    @debug "Computing next flip time."
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
    @debug "Next flip time is $smallestT."

    tTarget = 0

    if length(targets(T)) == 1
        tTarget = 1
    else
        tTarget = Nemo.PosInf()
    end

    if smallestT >= tTarget
        @debug "Flips happen after a target, so we need to check if we reached endgame or if a rebase needs to happen."
        if isinf(tTarget)
            @debug "Rebase needs to happen as we are sending a height to infinity."
            # set the heights to the first target heights and remove it from the list
            update_number_of_rebases!(T, 1)
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
        update_number_of_simultaneous_bergman_and_jensen_moves!(T, 1)
        @debug "Both Bergman and Jensen times are equal. Performing both moves."
    end

    changingBergmanMixedCells = [σ for σ in mixed_cells(T) if bergmanTimes[σ] == smallestT]
    changingJensenMixedCells = [σ for σ in mixed_cells(T) if jensenTimes[σ] == smallestT]

    @debug "Mixed cells that are changing: Bergman cells: $(changingBergmanMixedCells) and Jensen cells: $(changingJensenMixedCells)."

    # these two lists should have empty intersection, if not, then a perturbation was required
    @assert isempty(intersect(changingBergmanMixedCells, changingJensenMixedCells)) "A perturbation is required, but these are not implemented yet. Violating mixed cell(s): $(intersect(changingBergmanMixedCells, changingJensenMixedCells)). Bergman time(s): $(bergman_time.(Ref(T), intersect(changingBergmanMixedCells, changingJensenMixedCells))). Jensen time(s): $(jensen_time.(Ref(T), intersect(changingBergmanMixedCells, changingJensenMixedCells)))"

    newMixedCells = MixedCell[]

    if smallestTBergman == smallestT
        @debug "Performing Bergman move."
        update_number_of_bergman_moves!(T, 1)
        for σ in changingBergmanMixedCells
            @debug "Changing mixed cell: $σ."
            push!(newMixedCells, bergman_flip(T, σ, smallestTBergman)...)
            @debug "$(length(newMixedCells)) new mixed cells: $newMixedCells."
            remove_mixed_cell!(T, σ)
        end
    end
    l = length(newMixedCells)
    if smallestTJensen == smallestT
        @debug "Performing Jensen move."
        update_number_of_jensen_moves!(T, 1)
        for σ in changingJensenMixedCells
            @debug "Changing mixed cell: $σ."
            push!(newMixedCells, jensen_flip(T, σ, smallestTJensen)...)
            @debug "$(length(newMixedCells) - l) new mixed cells: $(newMixedCells[l+1:end])."
            remove_mixed_cell!(T, σ)
        end
    end

    # move the tracker to the heights achieved at the flip time
    @debug "Moving the tracker to the heights achieved at the flip time."
    add_heights!(T, smallestT*direction(T))

    # check consistency of new cells
    for σ in newMixedCells
        @assert is_transverse(σ) "$(σ) is not transverse"
        @assert are_support_heights_finite(T, σ) "$(σ) has invalid mixed height data"
        @assert is_bergman_consistent(T, σ) "Cells are not Bergman consistent"
    end

    @debug "New heights: $(dump_info(ambient_support(T)))."

    @debug "Merging $(length(newMixedCells)) new mixed cells."
    for σ in newMixedCells
        @debug "Merging mixed cell: $σ."
        merge_mixed_cell!(T, σ)
    end

    # update the cached jensen and bergman times for the mixed cells that didn't change
    update_cached_times!(T, newMixedCells, smallestT)

    @debug "Tracker move successfully."

    update_max_mixed_cells!(T, length(mixed_cells(T)))

    return true

end

function tropical_homotopy_continuation(T::Tracker)
    while move!(T) end
end

function stable_intersection(T::Tracker)::Vector{TropicalPoint}
    while move!(T) 
    @info "$(T.logger)"
    end
    @info "$(T.logger)"
    return [first(tropical_intersection_point_and_drift(T, σ)) for σ in mixed_cells(T)]
end