@doc raw"""
    starting_data(targetΔ::MixedSupport, M::RealisableMatroid)

Compute the starting data for the homotopy from the mixed support `targetΔ` and the realisable matroid `M`.

Returns a tuple (Δ, Δ', σ) with Δ a mixed support, Δ' the t arget, and σ the unique mixed cell in Δ.
"""
function starting_data(targetΔ::MixedSupport, M::RealisableMatroid)

    @vprintln :TropicalHomotopyContinuationStart "Computing starting data for tropical homotopies"
    @vprintln :TropicalHomotopyContinuationStart "Target mixed support: $(dump_info(targetΔ))"

    startingPolynomials = TropicalPolynomial[]
    TT = tropical_semiring()

    n = length(first(points(targetΔ)))
    R, x = polynomial_ring(TT, ["x$i" for i in 1:n])
    degrees = []
    for S in supports(targetΔ)

        @vprintln :TropicalHomotopyContinuationStart "Computing starting polynomial for point support $(S)"
        pts = points(S)
        ambientDim = length(first(pts))
        d = Int(ceil(max([norm(entries(p)) for p in pts]...)))
        push!(degrees, d)
        @vprintln :TropicalHomotopyContinuationStart "Degree d = $d"

        simplexVertices = Vector{Int}[]
        for j in 1:ambientDim
            @vprintln :TropicalHomotopyContinuationStart "Checking if $j-th simplex vertex is needed"
            if any([p[j] != 0 for p in pts])
                @vprintln :TropicalHomotopyContinuationStart "Adding $j-th simplex vertex $([j == i ? 1 : 0 for i in 1:ambientDim])"
                push!(simplexVertices, [i == j ? 1 : 0 for i in 1:ambientDim])
            else
                @vprintln :TropicalHomotopyContinuationStart "Simplex vertex $j is not needed"
            end
        end

        # work out whether to add 0 or not
        @vprintln :TropicalHomotopyContinuationStart "Checking if origin is needed"
        if any(norm(p) != d for p in entries.(pts))
            @vprintln :TropicalHomotopyContinuationStart "Adding origin"
            push!(simplexVertices, [0 for i in 1:ambientDim])
        else
            @vprintln :TropicalHomotopyContinuationStart "Origin is not needed"
        end

        # create a tropical linear polynomial using the simplexVertices as exponent vectors
        @vprintln :TropicalHomotopyContinuationStart "Creating tropical linear polynomial"
        f = 0
        for exponentVector in simplexVertices
            monomial = 1
            for i in 1:ambientDim
                monomial *= x[i]^exponentVector[i]
            end
            # make sure that the vertices of the simplex are lifted lower than everything inside targetΔ
            f += TT(QQ(min([targetΔ[p] for p in points(targetΔ)]...)) - QQ(rand(33:1000)))*monomial
        end
        @vprintln :TropicalHomotopyContinuationStart "Starting polynomial: $f"
        push!(startingPolynomials, f)

    end

    @vprintln :TropicalHomotopyContinuationStart "Finished computing all starting polynomials"
    @vprintln :TropicalHomotopyContinuationStart "Finding the tropical point corresponding to system $(startingPolynomials) with the linear equation matrix $(matrix(M))"
    w = find_tropical_point(startingPolynomials, M)
    @vprintln :TropicalHomotopyContinuationStart "Found tropical point: $w"
    C = chain_of_flats(M, QQ.(w))
    @vprintln :TropicalHomotopyContinuationStart "Chain of flats: $C"
    startingSupports = support.(startingPolynomials)
    # scale the entries of all the startingSupports by the degrees
    for (i, S) in enumerate(startingSupports)
        D = entries(S)
        newPoints = Point[]
        newHeights = QQFieldElem[]
        for p in keys(D)
            push!(newPoints, point(degrees[i] * p.entries))
            push!(newHeights, degrees[i] * D[p])
        end
        startingSupports[i] = support(newPoints, newHeights)
    end
    @vprintln :TropicalHomotopyContinuationStart "Active starting support: $(startingSupports)"
    # work out the active support for the mixed cell
    # by construction it will be monomials in startingSupports

    @vprintln :TropicalHomotopyContinuationStart "Computing active support for mixed cell"
    activeSupport = Support[]
    for S in startingSupports
        @assert length(minimum_monomials(S, QQ.(w))) == 2 "Not a valid starting support"
        push!(activeSupport, support(minimum_monomials(S,QQ.(w)), [0 for i in 1:length(minimum_monomials(S, QQ.(w)))]))
        @vprintln :TropicalHomotopyContinuationStart "Active support of mixed cell: $(activeSupport[end])"
    end

    # create a mixed cell from the active support
    @vprintln :TropicalHomotopyContinuationStart "Creating mixed cell from active support"
    σ = mixed_cell(mixed_support(tuple(activeSupport...)), C)
    @vprintln :TropicalHomotopyContinuationStart "Mixed cell: $σ"
    # merge targetΔ into startingSupports
    @vprintln :TropicalHomotopyContinuationStart "Merging target support into active starting support"
    for (i,S) in enumerate(startingSupports)
        @vprintln :TropicalHomotopyContinuationStart "Merging starting $(S)_$(i) with target support $(supports(targetΔ)[i])"
        startingSupports[i] = merge(S, supports(targetΔ)[i])
        @vprintln :TropicalHomotopyContinuationStart "Merged support: $(startingSupports[i])"
    end

    @vprintln :TropicalHomotopyContinuationStart "Creating mixed support from merged supports"
    Δ = mixed_support(tuple(startingSupports...))

    @vprintln :TropicalHomotopyContinuationStart "Mixed support: $(dump_info(Δ))"

    @vprintln :TropicalHomotopyContinuationStart "Finished computing starting data"

    if AbstractAlgebra.get_assertion_level(:TropicalHomotopyContinuationStart) > 0
        @assert is_subset(active_support(σ), Δ) "The active support of the mixed cell is not a subset of the mixed support."
        @assert is_transverse(σ) "The mixed cell is not transverse."
        @assert are_support_heights_finite(Δ, σ) "$(σ) has invalid mixed height data"
    end
    return Δ, σ
end

@doc raw"""
    find_tropical_point(f::TropicalPolynomial, M::RealisableMatroid)

Compute the tropical point in the tropical hypersurfaces of the elements of `polynomials` and the tropicalisation of the linear ideal associated to `M`.
"""
function find_tropical_point(polynomials::Vector{TropicalPolynomial}, M::RealisableMatroid)

    @assert all(total_degree(polynomial) == 1 for polynomial in polynomials) "The polynomials must be linear."

    R, t = rational_function_field(QQ, "t")
    S, _ = polynomial_ring(R, ngens(parent(first(polynomials))))

    ν = tropical_semiring_map(R, t)

    liftedPolynomials = random_lift.(Ref(ν), polynomials, Ref(S))

    rows = Vector{Oscar.elem_type(R)}[]
    for f in liftedPolynomials
        push!(rows, coeff.(Ref(f), gens(S)))
    end

    A = Oscar.matrix(R, rows)
    # for now, test whether this is working correctly so far
    linearIdealMatrix = Oscar.matrix(QQ, transpose(Oscar.nullspace(matrix(M))[2]))
    equationsMatrix = vcat(R.(matrix(M)), A)
    equationsMatrix = vcat(R.(linearIdealMatrix), A)

    # right hand side of the equations is 0 for the matroid parts
    b = zeros(R, nrows(equationsMatrix))
    for i in 1:nrows(linearIdealMatrix)
        b[i] = R(0)
    end
    for i in 1:length(liftedPolynomials)
        b[nrows(linearIdealMatrix) + i] = -constant_coefficient(liftedPolynomials[i])
    end

    try
        @vprintln :TropicalHomotopyContinuationStart "Solving system of equations with number of variables $(ncols(equationsMatrix)) and number of equations $(nrows(equationsMatrix))"
        w = Oscar.solve(equationsMatrix, b, side = :right)
        @vprintln :TropicalHomotopyContinuationStart "Found solution $(w)"
        return ν.(w)
    catch
        return nothing
    end

end

function random_lift(nu::TropicalSemiringMap, a::TropicalSemiringElem=tropical_semiring(nu)(rand(Int8)))

    aInt = ZZ(a; preserve_ordering=true)
    randomLift = rand(-999:999)*Oscar.uniformizer_field(nu)^aInt
    while iszero(randomLift)
        randomLift = rand(-999:999)*Oscar.uniformizer_field(nu)^aInt
    end
    return randomLift

end

@doc raw"""
    random_lift(nu::TropicalSemiringMap, f::TropicalPolynomial, parent)

Randomly lift each coefficient of `f` and return the result in the correct polynomial ring.
"""
function random_lift(nu::TropicalSemiringMap, f::TropicalPolynomial, parent)

    # random lift each coefficient of f, return result in correct polynomial_ring
    return map_coefficients(c -> random_lift(nu, c), f, parent=parent)

end
