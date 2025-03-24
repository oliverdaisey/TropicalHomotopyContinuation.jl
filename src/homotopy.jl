@doc raw"""
    move!(T::Tracker)

Move the heights of the tracker `T` to the next flip time, updating all mixed cells accordingly.

If we would reach the target before or at any flip time, do nothing.
"""
function move!(T::Tracker)
    
    # work out which mixed cell(s) have minimal Bergman / Jensen time
    smallestTBergman = minimum([bergman_time(T, σ) for σ in mixed_cells(T)])
    smallestTJensen = minimum([jensen_time(T, σ) for σ in mixed_cells(T)])
    smallestT = min(smallestTBergman, smallestTJensen) # the time at which we perform flips

    tTarget = 1

    if smallestT >= tTarget
        println("No moves possible before reaching the target mixed subdivision. End game reached.")
        return
    end

    changingBergmanMixedCells = [σ for σ in mixed_cells(T) if bergman_time(T, σ) == smallestT]
    changingJensenMixedCells = [σ for σ in mixed_cells(T) if jensen_time(T, σ) == smallestT]

    # these two lists should have empty intersection, if not, then a perturbation was required
    @assert isempty(intersect(changingBergmanMixedCells, changingJensenMixedCells)) "A perturbation is required, but these are not implemented yet. Violating mixed cell(s): $(intersect(changingBergmanMixedCells, changingJensenMixedCells)). Bergman time(s): $(bergman_time.(Ref(T), intersect(changingBergmanMixedCells, changingJensenMixedCells))). Jensen time(s): $(jensen_time.(Ref(T), intersect(changingBergmanMixedCells, changingJensenMixedCells)))"

    newMixedCells = MixedCell[]

    if smallestTBergman <= smallestTJensen
        for σ in changingBergmanMixedCells
            push!(newMixedCells, bergman_flip(T, σ)...)
            remove_mixed_cell!(T, σ)
        end
    end
    if smallestTJensen <= smallestTBergman
        for σ in changingJensenMixedCells
            push!(newMixedCells, jensen_flip(T, σ)...)
            remove_mixed_cell!(T, σ)
        end
    end

    # move the tracker to the heights achieved at the flip time
    add_heights!(T, min(smallestTBergman, smallestTJensen)*direction(T))

    for σ in newMixedCells
        merge_mixed_cell!(T, σ)
    end

end