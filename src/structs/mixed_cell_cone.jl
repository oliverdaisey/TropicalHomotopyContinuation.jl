@doc raw"""
    struct MixedCellConeFacet

A facet of a mixed cell cone defined by `circuit` that involves the extra point `point`.
"""
struct MixedCellConeFacet

    circuit::Dict{Point, Height}
    point::Point

end

@doc raw"""
    struct MixedCellCone

A mixed cell cone defined by `facets` over `ambientSupport`. These encode the permissible heights that give rise to the a mixed cell.
"""
struct MixedCellCone

    ambientSupport::MixedSupport
    facets::AbstractVector{MixedCellConeFacet}

end

@doc raw"""
    facets(C::MixedCellCone)

Return the facets of the mixed cell cone `C`.
"""
function facets(C::MixedCellCone)
    return C.facets
end

@doc raw"""
    extra_point(κ::MixedCellConeFacet)

Return the extra point outside the defining mixed cell candidate defining the mixed cell cone facet `κ`.
"""
function extra_point(κ::MixedCellConeFacet)
    return κ.point
end

function Base.show(io::IO, C::MixedCellCone)

    print(io, "Mixed cell cone with $(length(facets(C))) facets")
end

function (io::IO, F::MixedCellConeFacet)

    print(io, "Mixed cell cone facet corresponding to point $(extra_point(F))")
end

@doc raw"""
    mixed_cell_cone(facets::AbstractVector{MixedCellConeFacet})::MixedCellCone

Construct a mixed cell cone from `facets`.
"""
function mixed_cell_cone(ambientSupport::MixedSupport, facets::AbstractVector{MixedCellConeFacet})::MixedCellCone
    return MixedCellCone(ambientSupport, facets)
end

@doc raw"""
    mixed_cell_cone_facet(circuit::Dict{Point, height})::MixedCellConeFacet

Construct a mixed cell cone facet from `circuit` using the extra point `p`. These are the nontrivial entries of the defining linear functional of the facet.
"""
function mixed_cell_cone_facet(circuit::Dict{Point, Height}, p::Point)::MixedCellConeFacet
    return MixedCellConeFacet(circuit, p)
end

@doc raw"""
    circuit(κ::MixedCellConeFacet)

Return the circuit defining the mixed cell cone facet `κ`.
"""
function circuit(κ::MixedCellConeFacet)
    return κ.circuit
end

function Base.getindex(κ::MixedCellConeFacet, p::Point)
    # return the height of point p in the circuit defining the facet
    # return 0 if the mapping doesn't exist
    return get(circuit(κ), p, 0)
end

function mixed_cell_cone(δ::MixedSupport, ambientSupport::MixedSupport)::MixedCellCone

    @assert length(δ) == length(ambientSupport) "Mixed cell candidate and ambient support must have the same number of supports."
    @assert is_subset(δ, ambientSupport) "Mixed cell candidate must be a subset of the ambient support."

    # Compute the base Cayley matrix for delta and invert it once.
    # For each extra point p, the circuit of [B | c_p] is proportional to [-B^{-1} c_p; 1],
    # avoiding a full nullspace computation per point.
    baseMatrix = matrix(δ)
    B_qq = Oscar.matrix(QQ, baseMatrix)
    _, Binv = Oscar.is_invertible_with_inverse(B_qq)

    delta_pts = points(δ)
    numSupports = length(supports(δ))

    facets = MixedCellConeFacet[]

    for p in points(ambientSupport)
        if p in delta_pts
            continue
        end

        index = findfirst(x -> p in x, supports(ambientSupport))

        # Build the Cayley column for p: [entries(p); e_index]
        col = vcat(entries(p), [i == index ? 1 : 0 for i in 1:numSupports])
        c_qq = Oscar.matrix(QQ, length(col), 1, col)

        # Circuit entries for delta points are -B^{-1} * c, and 1 for p
        x = Binv * (-c_qq)

        circuit = Dict{Point, Height}()
        for (i, point) in enumerate(delta_pts)
            circuit[point] = x[i, 1]
        end
        circuit[p] = QQ(1)

        # Choose sign so that the entry corresponding to p is negative
        if circuit[p] > 0
            for pt in keys(circuit)
                circuit[pt] = -circuit[pt]
            end
        end

        push!(facets, mixed_cell_cone_facet(circuit, p))

    end

    return mixed_cell_cone(ambientSupport, facets)
end

function Base.show(io::IO, C::MixedCellConeFacet)

    # print each support with their active points
    println("Mixed cell cone facet with extra point $(extra_point(C))")
end

@doc raw"""
    mixed_cell_cone(candidate::MixedSupport, ambientSupport::MixedSupport)

Compute the mixed cell cone of a mixed cell candidate `σ` with ambient support `ambientSupport`.
"""
function mixed_cell_cone(σ::MixedCell, ambientSupport::MixedSupport)::MixedCellCone

    return mixed_cell_cone(active_support(σ), ambientSupport)
end

@doc raw"""
    in(Δ::MixedSupport, C::MixedCellCone)

Tests whether the heights in `Δ` are in the mixed cell cone `C`.
"""
function Base.in(Δ::MixedSupport, C::MixedCellCone)
    pts = points(Δ)

    dotProducts = []
    for κ in facets(C)
        push!(dotProducts, sum([QQ(circuit(κ)[p]) * QQ(Δ[p]) for p in pts if p in keys(circuit(κ))]))
    end

    if all(dotProducts .<= 0)
        return true
    else
        return false
    end
end

@doc raw"""
    ambient_support(C::MixedCellCone)

Return the ambient support defining the mixed cell cone `C`.
"""
function ambient_support(C::MixedCellCone)
    return C.ambientSupport
end

@doc raw"""
    points(C::MixedCellCone)

Return the ambient points that define the support of the mixed cell cone `C`.
"""
function points(C::MixedCellCone)
    return points(ambient_support(C))
end

@doc raw"""
    Base.convert(::Type{Polyhedron}, C::MixedCellCone)

Convert a mixed cell cone `C` to a polymake polyhedron.
"""
function Base.convert(::Type{Polyhedron}, C::MixedCellCone)
    pts = points(C)
    A = Vector{QQFieldElem}[]
    b = QQFieldElem[]
    for κ in facets(C)
        push!(A, [p in keys(circuit(κ)) ? circuit(κ)[p] : 0 for p in pts])
        push!(b, 0)
    end

    # convert A to polymake compatible format

    A = Oscar.matrix(QQ, A)

    return polyhedron(A, b)
end

function dot(Δ::MixedSupport, κ::MixedCellConeFacet)
    return sum([QQ(circuit(κ)[p]) * QQ(Δ[p]) for p in keys(circuit(κ)) if !isinf(Δ[p])])
end

function dot(κ::MixedCellConeFacet, Δ::MixedSupport)
    return dot(Δ, κ)
end

function dot(Δ::MixedSupport, Ε::MixedSupport)
    pts = points(Δ)
    return sum([Δ[p] * Ε[p] for p in pts])
end

function Base.in(κ::MixedCellConeFacet, Δ::MixedSupport)
    return dot(Δ, κ) == 0
end
