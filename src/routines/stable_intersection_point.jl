
"""
    stable_intersection_point(s::MixedCell)

Returns the tropical intersection point dual to the mixed cell `s`.
"""
function stable_intersection_point(s::MixedCell)

    # we will construct the affine equations that determines the intersection point
    exponentDifferences = Vector{Vector{Int}}()
    coefficientDifferences = Vector{Oscar.TropicalSemiringElem}()
    for si in dual_cells(s)
        activeSupport = active_support(si)
        dualWeight = dual_weight(si)
        # only consider exponents with tropically nonzero coefficients
        nonTrivialPoints = filter(x -> !iszero(dualWeight[x]), activeSupport)
        # iterate through all 2-element subsets of tropicallyNonZeroIndices
        for subset in subsets(nonTrivialPoints, 2)
            k, l = subset
            # coefficientDifferences of the new equation
            push!(exponentDifferences, k-l)
            # constant of the new equation
            push!(coefficientDifferences, dualWeight[l] / dualWeight[k])
        end
    end

    M = matrix(QQ, exponentDifferences)
    b = QQ.(coefficientDifferences; preserve_ordering=true)
    @assert rank(M) == ncols(M) - tropical_lineality_dim(dual_cells(s)) "The intersection is not transverse modulo lineality"

    linealityGens = tropical_lineality_space(s)
    M = vcat(M, linealityGens)
    b = vcat(b, zeros(QQ, nrows(linealityGens)))

    solution = solve(M, b, side=:right)
    return solution

end