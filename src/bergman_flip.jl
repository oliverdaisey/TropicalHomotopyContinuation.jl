@doc raw"""
    compute_bergman_time(T::Tracker)

Compute the Bergman time of the mixed cell `σ` with tracker `T`. This is the smallest time t at which `w + t * u` induces a coarser chain of flats than the current chain of flats, where `w` and `u` are the tropical intersection point and tropical drift of `T` respectively.

Caches the answer inside `T`.
"""
function compute_bergman_time(T::Tracker, σ::MixedCell; disable_cache::Bool = false)

    chainOfFlats = chain_of_flats(σ)
    w, u = tropical_intersection_point_and_drift(T, σ)
    inequalities, equalities = cone(chainOfFlats)

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationBergman)>0
        # check that we are inside the cone
        if length(equalities) > 0
            @assert all([sum(w .* v) == 0 for v in equalities]) "The intersection point is not in the cone (equality violated)"
            @assert all([sum(w .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone (inequality violated) intersection point: $(w) chain of flats: $(chainOfFlats)"
        else
            @assert all([sum(w .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone"
        end
    end

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationBergman)>1
        @assert u in cone_from_equations(linear_equation_matrix(linear_span(C))) "The drift is not in the cone"
    end

    timesOfIntersection = [sum(v.*u) != 0 ? -sum(v.*w) / sum(v.*u) : Nemo.PosInf() for v in inequalities]
    # delete all times that are less than 0
    timesOfIntersection = [t for t in timesOfIntersection if t > 0]
    if timesOfIntersection == []
        return Nemo.PosInf()
    end

    # check that we are still inside the cone
    t = minimum(timesOfIntersection)
    if !isinf(t)
        if length(equalities) > 0
            @assert all([sum((w + t*u) .* v) == 0 for v in equalities]) "The intersection point is not in the cone (equality violated)"
        end
        for v in inequalities
            if sum((w + t*u) .* v) > 0
                @assert false "Inequality $(v) is violated at $(w + t*u) corresponding to t = $(t). The time of intersection for this facet is equal to $(sum(v.*u) != 0 ? -sum(v.*w) / sum(v.*u) : Nemo.PosInf())"
            end
        end
        @assert all([sum((w + t*u) .* v) <= 0 for v in inequalities]) "The intersection point is not in the cone (inequality violated) intersection point: $(w) chain of flats: $(chainOfFlats)"
    end

    bergmanTime = minimum(timesOfIntersection)

    # cache the answer
    if !disable_cache
        update_bergman_time!(T, σ, bergmanTime)
    end

    return bergmanTime
end

@doc raw"""
    bergman_flip(T::Tracker)

Compute the Bergman flip of the mixed cell `σ` with tracker `T`.
"""
function bergman_flip(T::Tracker, σ::MixedCell, tBergman::Height)

    # Compute the tropical intersection point and drift
    w, u = tropical_intersection_point_and_drift(T, σ)

    # Compute the unrefined chain at the point of the Bergman flip
    startingChain = chain_of_flats(σ)
    unrefinedChain = chain_of_flats(matroid(startingChain), w + tBergman * u)
    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationBergman)>0
        # check that the unrefined chain is indeed non-maximal
        # if violated, then input tBergman does not match computed w and u
        @assert length(unrefinedChain) < length(startingChain) "The unrefined chain is maximal"
    end

    # If unrefined chain has colength greater 1, we require a perturbation
    if length(unrefinedChain) < + 1 != length(startingChain)
        return MixedCell[] # empty return signal need for perturbation
    end

    # Otherwise, compute all maximal refinements
    refinedChains = maximal_refinements(unrefinedChain)

    # Optimisation: remove original chain
    # refinedChains = [chain for chain in refinedChains if chain != C]

    # Next, compute the jensen trail, i.e., we consider the stars of the tropical hypersurfaces
    # around the intersection point, which only depends on the active support, and construct the
    # path their intersection takes with varying t
    jensenStarEquations = Vector{QQFieldElem}[]
    for Si in supports(active_support(σ))
        pts = points(Si)
        p1 = first(pts) # TODO: can we use Iterators.peel here instead of isequal below?
        for p in pts
            if !is_equal(p1, p)
                push!(jensenStarEquations, convert(Vector{QQFieldElem}, p1 - p))
            end
        end
    end
    jensenStarLinearEquationMatrix = Oscar.matrix(QQ, jensenStarEquations)
    jensenStarAffineEquationsRHS = jensenStarLinearEquationMatrix * (w + (tBergman+1) * u)

    # ###
    # # Old test
    # ###
    # jensenTrail = convex_hull([w + tBergman * u], [u], transpose(kernel(jensenStarLinearEquationMatrix, side=:right))) # only needed for old test
    # flippedChains = ChainOfFlats[]
    # for refinedChain in refinedChains
    #     # construct the rays of the bergmanCone
    #     bergmanConeRays = [ indicator_vector(Fj) for Fj in full_flats(refinedChain) if !isempty(Fj) ]
    #     bergmanConeRayMatrix = Oscar.matrix(QQ, bergmanConeRays)

    #     # quick test: is the intersection of the bergmanCone and jensenStar transverse,
    #     # i.e., does a bergmanConeRay lie in jensenStar?
    #     if Oscar.rank(bergmanConeRayMatrix) != Oscar.rank(jensenStarLinearEquationMatrix*transpose(bergmanConeRayMatrix))
    #         continue
    #     end

    #     # comprehensive test: does the jensenTrail intersect bergmanCone in a line?
    #     chainOfFlatsCone = positive_hull(bergmanConeRayMatrix, ones_matrix(QQ, 1,ncols(bergmanConeRayMatrix)))
    #     if Oscar.dim(intersect(jensenTrail, polyhedron(chainOfFlatsCone))) != 1
    #         continue
    #     end

    #     push!(flippedChains, refinedChain)
    # end

    ###
    # New test
    ###
    # newFlippedChains = ChainOfFlats[]
    flippedChains = ChainOfFlats[]
    for refinedChain in refinedChains
        # construct the linear equations of the span of the Bergman cone
        bergmanConeRays = [ indicator_vector(Fj) for Fj in full_flats(refinedChain) if !isempty(Fj) ]
        bergmanConeRayColumnMatrix = transpose(Oscar.matrix(QQ, bergmanConeRays))
        bergmanConeSpanLinearEquations = kernel(bergmanConeRayColumnMatrix)

        # assemble and solve the combined affine linear system from the Bergman linear equations
        # and the jensen affine linear equations
        affineLinearEquationsLHS = vcat(bergmanConeSpanLinearEquations, jensenStarLinearEquationMatrix)
        affineLinearEquationsRHS = vcat(zeros(QQ, nrows(bergmanConeSpanLinearEquations)), jensenStarAffineEquationsRHS)
        canSolve, solution, kernelGenerators = Oscar.can_solve_with_solution_and_kernel(affineLinearEquationsLHS, affineLinearEquationsRHS; side=:right)

        # check whether the combined affine linear system has exactly one solution.
        # if solution set is empty or positive-dimensional, then Jensen star and Bergman cone
        # do not intersect transversally
        if !canSolve || ncols(kernelGenerators)>0
            continue
        end

        # check whether solution induces refinedChain from unrefinedChain
        if !(induces_refinement(refinedChain, unrefinedChain, solution))
            continue
        end

        push!(flippedChains, refinedChain)
        # push!(newFlippedChains, refinedChain)
    end

    # println("=======================================================")
    # println("=======================================================")
    # println("=======================================================")

    # if !issetequal(newFlippedChains, flippedChains)
    #     println("Unrefined chain: $(unrefinedChain)")
    #     println("Old flipped chains: $(flippedChains)")
    #     println("New flipped chains: $(newFlippedChains)")

    #     println("w: $(w)")
    #     println("u: $(u)")
    #     println("tBergman: $(tBergman)")

    #     println("jensenStarLinearEquationMatrix: $(jensenStarLinearEquationMatrix)")

    #     ###
    #     # Debugging first old flipped chain that is not in the new flipped chains
    #     ###
    #     refinedChain = first(flippedChains)
    #     println("taking first old refined chain: $(refinedChain)")

    #     bergmanConeRays = [ indicator_vector(Fj) for Fj in full_flats(refinedChain) if !isempty(Fj) ]
    #     bergmanConeRayColumnMatrix = transpose(Oscar.matrix(QQ, bergmanConeRays))
    #     bergmanConeSpanLinearEquations = kernel(bergmanConeRayColumnMatrix)
    #     println("linear equations of the span of the Bergman cone:")
    #     display(bergmanConeSpanLinearEquations)

    #     affineLinearEquationsLHS = vcat(bergmanConeSpanLinearEquations, jensenStarLinearEquationMatrix)
    #     affineLinearEquationsRHS = vcat(zeros(QQ, nrows(bergmanConeSpanLinearEquations)), jensenStarAffineEquationsRHS)
    #     println("bergmanConeSpanLinearEquations: $(bergmanConeSpanLinearEquations)")
    #     display(bergmanConeSpanLinearEquations)
    #     println("jensenStarLinearEquationMatrix: $(jensenStarLinearEquationMatrix)")
    #     display(jensenStarLinearEquationMatrix)
    #     println("jensenStarAffineEquationsRHS: $(jensenStarAffineEquationsRHS)")
    #     println("affine linear equations LHS: $(affineLinearEquationsLHS)")
    #     display(affineLinearEquationsLHS)
    #     println("affine linear equations RHS: $(affineLinearEquationsRHS)")

    #     canSolve, solution, kernelGenerators = Oscar.can_solve_with_solution_and_kernel(affineLinearEquationsLHS, affineLinearEquationsRHS; side=:right)
    #     println("can solve: $(canSolve)")
    #     println("solution: $(solution)")
    #     println("kernel generators: $(kernelGenerators)")

    #     bergmanCone = positive_hull(transpose(bergmanConeRayColumnMatrix), ones_matrix(QQ, 1,nrows(bergmanConeRayColumnMatrix)))
    #     println("bergmanConeRayColumnMatrix: $(bergmanConeRayColumnMatrix)")
    #     println("solution in bergmanCone: $(solution in bergmanCone)")

    #     println("breaking direction: $(breaking_direction(refinedChain, unrefinedChain))")
    #     println("bergmanConeSpanLinearEquations: $(bergmanConeSpanLinearEquations)")
    #     println("solution: $(solution)")

    #     # ###
    #     # # Debugging last new flipped chain not in old flipped chains
    #     # ###
    #     # refinedChain = last(newFlippedChains)
    #     # println("taking last new refined chain: $(refinedChain)")

    #     # bergmanConeRays = [ indicator_vector(Fj) for Fj in full_flats(refinedChain) if !isempty(Fj) ]
    #     # bergmanConeRayMatrix = Oscar.matrix(QQ, bergmanConeRays)
    #     # println("bergmanConeRayMatrix: $(bergmanConeRayMatrix)")

    #     # println("Oscar.rank(bergmanConeRayMatrix): $(Oscar.rank(bergmanConeRayMatrix))")
    #     # println("Oscar.rank(jensenStarLinearEquationMatrix*transpose(bergmanConeRayMatrix)): $(Oscar.rank(jensenStarLinearEquationMatrix*transpose(bergmanConeRayMatrix)))")

    #     # chainOfFlatsCone = positive_hull(bergmanConeRayMatrix, ones_matrix(QQ, 1,ncols(bergmanConeRayMatrix)))
    #     # println("Oscar.dim(intersect(jensenTrail, polyhedron(chainOfFlatsCone))): $(Oscar.dim(intersect(jensenTrail, polyhedron(chainOfFlatsCone))))")

    #     @assert false "The new flipped chains and the old flipped chains should be the same"
    # end


    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationBergman)>0
        # check that the mixed cell data is valid
        newMixedCells = mixed_cell.(Ref(active_support(σ)), flippedChains)
        for σ in newMixedCells
            # check that the matrix coming from σ is invertible
            @assert is_transverse(σ) "$(σ) is not transverse"
            @assert are_support_heights_finite(T, σ) "$(σ) has invalid mixed height data"
        end
        @assert length(newMixedCells) > 0 "No new mixed cells during a Bergman flip"
    end

    return mixed_cell.(Ref(active_support(σ)), flippedChains)

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

# orthogonal projection away from subspace spanned by columns of kernelMatrix
function projection_matrix(kernelMatrix)
    _, invKernelMatrix = Oscar.is_invertible_with_inverse(transpose(kernelMatrix) * kernelMatrix)
    return kernelMatrix * invKernelMatrix * transpose(kernelMatrix)
end

@doc raw"""
    oblique_projection_matrix(A::Matrix, B::Matrix)

Compute the oblique projection matrix from the subspace spanned by the columns of `A` onto the subspace spanned by the columns of `B`.
"""
function oblique_projection_matrix(A, B)
    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationBergman)>0
        @assert Oscar.rank(B) == ncols(B) "Matrix B is not full rank"
        @assert Oscar.rank(A) == ncols(A) "Matrix A is not full rank"
    end
    _, invKernelMatrix = Oscar.is_invertible_with_inverse(transpose(B) * A)
    return A*(invKernelMatrix)*transpose(B)
end
