@doc raw"""
    bergman_time(T::Tracker)

Compute the Bergman time of the mixed cell `σ` with tracker `T`. This is the smallest time t at which `w + t * u` induces a coarser chain of flats than the current chain of flats, where `w` and `u` are the tropical intersection point and tropical drift of `T` respectively.
"""
function bergman_time(T::Tracker, σ::MixedCell)

    chainOfFlats = chain_of_flats(σ)
    w, u = tropical_intersection_point_and_drift(T, σ)
    C = cone(chainOfFlats)

    @assert w in C "The intersection point is not in the cone"
    @assert u in cone_from_equations(linear_equation_matrix(linear_span(C))) "The drift is not in the cone"

    inequalities = linear_inequality_matrix(Oscar.facets(C))
    inequalities = Vector{QQFieldElem}[inequalities[i, :] for i in 1:nrows(inequalities)]
    return minimum([dot(v, u) != 0 ? -dot(v, w) / dot(v, u) : PosInf for v in inequalities])
end

@doc raw"""
    bergman_flip(T::Tracker)

Compute the Bergman flip of the mixed cell `σ` with tracker `T`.
"""
function bergman_flip(T::Tracker, σ::MixedCell)

    tBergman = bergman_time(T, σ)
    C = chain_of_flats(σ)
    Δ = ambient_support(T)
    w, u = tropical_intersection_point_and_drift(T, σ)
    
    refinedChains = maximal_refinements(chain_of_flats(matroid(C), w + tBergman * u))
    # make sure not to include the original chain
    refinedChains = [chain for chain in refinedChains if chain != C]

    rows = Vector{QQFieldElem}[]
    allowedChains = ChainOfFlats[]
    for chain in refinedChains
        # add rows for each indicator vector of chain
        push!(rows, indicator_vector.(flats(chain))...)

        for S in supports(Δ)
            pts = points(S)
            p1 = first(pts)
            localRows = Vector{QQFieldElem}[]
            for p in pts
                if !isequal(p1, p)
                    push!(localRows, convert(Vector{QQFieldElem}, p1 - p))
                end
            end
            # compute kernel of matrix given by rows
            _, kernel = Oscar.nullspace(Oscar.matrix(QQ, localRows))
            # push the columns of `kernel` into `rows`

            # first convert kernel to list of vectors
            kernel = Vector{QQFieldElem}[kernel[:, i] for i in 1:ncols(kernel)]
            # then push each vector into `rows`
            push!(rows, kernel...)
        end

        # check that `rows` has rank nonzero
        if Oscar.rank(Oscar.matrix(QQ, rows)) == 0
            continue
        end

        # check that dot product of sum of indicator vectors with drift is non zero
        if dot(u, sum(indicator_vector.(flats(chain)))) <= 0
            continue
        end

        push!(allowedChains, chain)
        
    end

    return mixed_cell.(Ref(active_support(σ)), allowedChains)

end

@doc raw"""
    bergman_move!(T::Tracker)

Perform a Bergman move on the tracker `T`. This updates the mixed heights and the mixed cells tracked by `T` to the point where a maximal cone of the fine structure is breached, assuming no maximal tropical polyhedra corresponding to the hypersurface dual cells are breached.

We can perform a Bergman move when `bergman_time(T)` is less than `jensen_time(T)`.
"""
function bergman_move!(T::Tracker)
    
    # work out which mixed cell(s) have minimal Bergman time
    smallestTBergman = minimum([bergman_time(T, σ) for σ in mixed_cells(T)])
    changingMixedCells = [σ for σ in mixed_cells(T) if bergman_time(T, σ) == smallestTBergman]

    newMixedCells = MixedCell[]
    for σ in changingMixedCells
        push!(newMixedCells, bergman_flip(T, σ)...)
        remove_mixed_cell!(T, σ)
    end

    # move the tracker to the heights achieved at the bergman time
    add_heights!(T, smallestTBergman*direction(T))

    for σ in newMixedCells
        merge_mixed_cell!(T, σ)
    end

end