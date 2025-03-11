using Oscar

function jensen_time(T::Tracker, σ::MixedCell)::Union{QQFieldElem,Nothing}
    
    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    C = mixed_cell_cone(δ, Δ)
    
    @assert δ in C "The mixed cell being tracked is not in the mixed cell cone."

    v = direction(T)

    return minimum([dot(v, κ) != 0 ? -dot(v, δ) / dot(v, κ) : PosInf for κ in facets(C)])
    
end


function jensen_flip(T::Tracker, σ::MixedCell)
    # work out which facet inequality gets broken
    tJensen = jensen_time(T, σ)

    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    v = direction(T)

    C = mixed_cell_cone(δ, Δ)

    @assert sum([dot(v, κ) * tJensen == -dot(v, κ) for κ in facets(C)]) == 1 "The mixed cell being tracked does not breach its mixed cell cone in exactly one facet."

    κ = findfirst(κ -> dot(v, κ) * tJensen == -dot(v, κ), facets(C))
    p = extra_point(κ)

    newMixedCells = MixedCell[]
    
    for q in points(δ)
        if κ[q] > 0
            # return mixed cell with q in support and p not in support
            push!(newMixedCells, swap(σ, p, q))
        end
    end

    for σ in newMixedCells
        push!(T.mixedCells, σ)
    end

    # change the heights by jensen_time*direction
    add_heights!(T, direction(T) * tJensen)
end

function jensen_move!(T::Tracker)
    
    # work out which mixed cell(s) have minimal Bergman time
    smallestTJensen = minimum([jensen_time(T, σ) for σ in mixed_cells(T)])
    changingMixedCells = [σ for σ in mixed_cells(T) if jensen_time(T, σ) == smallestTJensen]

    newMixedCells = MixedCell[]
    for σ in changingMixedCells
        push!(newMixedCells, jensen_flip(T, σ)...)
        remove_mixed_cell!(T, σ)
    end

    # move the tracker to the heights achieved at the bergman time
    add_heights!(T, smallestTJensen*direction(T))

    for σ in newMixedCells
        merge_mixed_cell!(T, σ)
    end

end