@doc raw"""
    starting_data(targetΔ::MixedSupport, M::RealisableMatroid)

Compute the starting data for the homotopy from the mixed support `targetΔ` and the realisable matroid `M`.

Returns a tuple (Δ, σ) with Δ a mixed support and σ the unique mixed cell in Δ.
"""
function starting_data(targetΔ::MixedSupport, M::RealisableMatroid)

    startingPolynomials = TropicalPolynomial[]
    TT = tropical_semiring()

    n = length(first(points(targetΔ)))
    R, x = polynomial_ring(TT, ["x$i" for i in 1:n])
    for S in supports(targetΔ)

        pts = points(S)
        ambientDim = length(first(pts))
        d = Int(ceil(max([norm(entries(p)) for p in pts]...)))

        simplexVertices = Vector{Int}[]
        for j in 1:ambientDim
            if any([p[j] != 0 for p in pts])
                push!(simplexVertices, [i == j ? 1 : 0 for i in 1:ambientDim])
            end
        end

        # work out whether to add 0 or not
        if any(norm(p) != d for p in entries.(pts))
            push!(simplexVertices, [0 for i in 1:ambientDim])
        end

        # create a tropical linear polynomial using the simplexVertices as exponent vectors
        f = 0
        for exponentVector in simplexVertices
            monomial = 1
            for i in 1:ambientDim
                monomial *= x[i]^exponentVector[i]
            end
            # make sure that the vertices of the simplex are lifted lower than everything inside targetΔ
            f += TT(min([targetΔ[p] for p in points(targetΔ)]...) - QQ(rand(1:1000)))*monomial
        end
        push!(startingPolynomials, f)

    end

    C = chain_of_flats(M, QQ.(find_tropical_point(startingPolynomials, M)))

    startingSupports = support.(startingPolynomials)

    # merge targetΔ into startingSupports
    for (i,S) in enumerate(startingSupports)
        startingSupports[i] = merge(S, supports(targetΔ)[i])
    end

    Δ = mixed_support(tuple(startingSupports...))
    return Δ, mixed_cell(Δ, C)
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
    equationsMatrix = vcat(R.(matrix(M)), A)

    # right hand side of the equations is 0 for the matroid parts
    b = zeros(R, nrows(equationsMatrix))
    for i in 1:nrows(matrix(M))
        b[i] = R(0)
    end
    for i in 1:length(liftedPolynomials)
        b[nrows(matrix(M)) + i] = -constant_coefficient(liftedPolynomials[i])
    end

    try
        w = Oscar.solve(equationsMatrix, b, side = :right)
        return ν.(w)
    catch
        return nothing
    end

end

function random_lift(nu::TropicalSemiringMap, a::Union{TropicalSemiringElem, Nothing}=nothing)
    functionField = Oscar.valued_field(nu)
    coefficientField = base_ring(functionField)

    if isnothing(a)
        a = ZZ(rand(Int8))
    else
        a = ZZ(a; preserve_ordering=true)
    end
    randomLift = functionField(uniformizer(nu))^a


    if coefficientField==QQ
        return rand(-999:999)*randomLift
    else
        return rand(K)*randomLift
    end
end

@doc raw"""
    random_lift(nu::TropicalSemiringMap, f::TropicalPolynomial, parent)

Randomly lift each coefficient of `f` and return the result in the correct polynomial ring.
"""
function random_lift(nu::TropicalSemiringMap, f::TropicalPolynomial, parent)

    # random lift each coefficient of f, return result in correct polynomial_ring
    return map_coefficients(c -> random_lift(nu, c), f, parent=parent)

end