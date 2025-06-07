export stable_intersection

"""
    stable_intersection(S::MixedSupport, M::RealisableMatroid)

Computes the stable intersection given the data of a mixed support `S` and a realisable matroid `M`.
"""
function stable_intersection(S::MixedSupport, M::RealisableMatroid)
    Δ, σ = starting_data(S, M)
    T = tracker(Δ, S, [σ], path=:coefficient_wise)

    return TropicalHomotopyContinuation.stable_intersection(T)
end

"""
    compute_mixed_support(v::Vector{TropicalPolynomial})

Computes the mixed support associated to a vector of polynomials `v`.
"""
function compute_mixed_support(v::Vector)
    R = parent(first(v))
    x = gens(R)

    supports = Support[]
    for f in v
        push!(supports, support(Point.(collect(exponents(f))), collect(coefficients(f))))
    end

    return mixed_support(supports)
end

import Oscar.stable_intersection
stable_intersection(v::Vector, M::RealisableMatroid) = stable_intersection(compute_mixed_support(v), M)