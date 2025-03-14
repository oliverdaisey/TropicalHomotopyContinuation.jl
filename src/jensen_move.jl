using Oscar

function jensen_time(T::Tracker, σ::MixedCell)::Union{QQFieldElem,PosInf}
    
    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    C = mixed_cell_cone(δ, Δ)

    @assert Δ in C "The mixed cell being tracked is not in the mixed cell cone."

    v = direction(T)

    timesOfIntersection = [dot(v, κ) != 0 ? -dot(κ, Δ) / dot(v, κ) : Nemo.PosInf() for κ in facets(C)]
    # delete all times that are less than 0
    timesOfIntersection = [t for t in timesOfIntersection if t > 0]
    if timesOfIntersection == []
        return Nemo.PosInf()
    end
    display(timesOfIntersection)
    return minimum(timesOfIntersection)
    
end


function jensen_flip(T::Tracker, σ::MixedCell)
    # work out which facet inequality gets broken
    tJensen = jensen_time(T, σ)

    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    v = direction(T)

    C = mixed_cell_cone(δ, Δ)

    @assert sum([dot(v, κ) * tJensen == -dot(Δ, κ) for κ in facets(C)]) == 1 "The mixed cell being tracked does not breach its mixed cell cone in exactly one facet."

    κ = facets(C)[findfirst(κ -> dot(v, κ) * tJensen == -dot(Δ, κ), facets(C))]
    p = extra_point(κ)
    changingSupport = supports(Δ)[findfirst(p in points(S) for S in supports(Δ))]
    println("changingSupport = ", changingSupport)
    # work out which mixed cell support this corresponds to (which active support)
    changingDualCell = supports(σ)[findfirst(is_subset(S, changingSupport) for S in supports(σ))]
    println("changingDualCell = ", changingDualCell)

    @assert length(changingDualCell) == 2 "The changing dual cells is not minimal."
    println("extra point = ", p)
    newMixedCells = MixedCell[]

    for q in points(changingDualCell)
        if κ[q] > 0
            # return mixed cell with p in support and q not in support
            push!(newMixedCells, swap(σ, q, p))
        end
    end

    return newMixedCells
end

function jensen_move!(T::Tracker)
    
    # work out which mixed cell(s) have minimal Bergman time
    smallestTJensen = minimum([jensen_time(T, σ) for σ in mixed_cells(T)])
    changingMixedCells = [σ for σ in mixed_cells(T) if jensen_time(T, σ) == smallestTJensen]

    println("Changing $(length(changingMixedCells)) mixed cells")

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