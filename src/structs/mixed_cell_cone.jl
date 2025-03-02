using Oscar

"""
    struct MixedCellConeFacet

A facet of a mixed cell cone defined by `circuit` that involves the extra point `point`.
"""
struct MixedCellConeFacet

    circuit::Dict{Point, Weight}
    point::Point

end

"""
    struct MixedCellCone

A mixed cell cone defined by `facets` over `ambientSupport`. These encode the permissible heights that give rise to the a mixed cell.
"""
struct MixedCellCone

    ambientSupport::MixedSupport
    facets::AbstractVector{MixedCellConeFacet}

end

"""
    facets(C::MixedCellCone)

Return the facets of the mixed cell cone `C`.
"""
function facets(C::MixedCellCone)
    return C.facets
end

"""
    extra_point(κ::MixedCellConeFacet)

Return the extra point outside the defining mixed cell candidate defining the mixed cell cone facet `κ`.
"""
function extra_point(κ::MixedCellConeFacet)
    return κ.point
end

function Base.show(io::IO, C::MixedCellCone)

    print(io, "Mixed cell cone with $(length(facets(C))) facets")
end

function Base.show(io::IO, F::MixedCellConeFacet)

    print(io, "Mixed cell cone facet corresponding to point $(extra_point(F))")
end

"""
    mixed_cell_cone(facets::AbstractVector{MixedCellConeFacet})::MixedCellCone

Construct a mixed cell cone from `facets`.
"""
function mixed_cell_cone(ambientSupport::MixedSupport, facets::AbstractVector{MixedCellConeFacet})::MixedCellCone
    return MixedCellCone(ambientSupport, facets)
end

"""
    mixed_cell_cone_facet(circuit::Dict{Point, Weight})::MixedCellConeFacet

Construct a mixed cell cone facet from `circuit` using the extra point `p`. These are the nontrivial entries of the defining linear functional of the facet.
"""
function mixed_cell_cone_facet(circuit::Dict{Point, Weight}, p::Point)::MixedCellConeFacet
    return MixedCellConeFacet(circuit, p)
end

"""
    circuit(κ::MixedCellConeFacet)

Return the circuit defining the mixed cell cone facet `κ`.
"""
function circuit(κ::MixedCellConeFacet)
    return κ.circuit
end

"""
    mixed_cell_cone(candidate::MixedSupport, ambientSupport::MixedSupport)

Compute the mixed cell cone of a mixed cell candidate `candidate` with ambient support `ambientSupport`.
"""
function mixed_cell_cone(candidate::MixedSupport, ambientSupport::MixedSupport)::MixedCellCone

    @assert length(candidate) == length(ambientSupport) "Mixed vell candidate and ambient support must have the same number of supports."

    cayleyEmbedding = cayley_embedding(ambientSupport)

    facets = MixedCellConeFacet[]
    
    for p in points(ambientSupport)
        if p in points(candidate)
            continue
        end
        # for all points p not in ambient support, get submatrix of cayleyEmbedding indexed by candidate and p
        index = findfirst(x -> p in x, supports(ambientSupport))

        # take the support with index `index` and augment it with p
        oldSupport = supports(candidate)[index]
        newSupport = support(oldSupport, p)
        newMixedSupport = mixed_support(candidate, oldSupport, newSupport)
        submatrix = cayleyEmbedding[newMixedSupport]
        nontrivialEntries = Matrix(nullspace(Oscar.matrix(QQ, submatrix))[2])

        # choose sign so that the entry corresponding to p is negative
        if nontrivialEntries[findfirst(x -> x == p, points(newMixedSupport))] > 0
            nontrivialEntries = -nontrivialEntries
        end
        
        # circuit has enties all zero except for nontrivialEntries
        circuit = Dict{Point, Weight}()
        for (point, value) in zip(points(candidate), nontrivialEntries)
            circuit[point] = weight(value)
        end

        push!(facets, mixed_cell_cone_facet(circuit, p))

    end

    return mixed_cell_cone(ambientSupport, facets)
end

"""
    in(Δ::MixedSupport, C::MixedCellCone)

Tests whether the weights in `Δ` are in the mixed cell cone `C`.
"""
function Base.in(Δ::MixedSupport, C::MixedCellCone)
    pts = points(Δ)

    dotProducts = []
    for κ in facets(C)
        push!(dotProducts, sum([circuit(κ)[p] * Δ[p] for p in pts if p in keys(circuit(κ))]))
    end

    return all(dotProducts .< 0)
end

"""
    ambient_support(C::MixedCellCone)

Return the ambient support defining the mixed cell cone `C`.
"""
function ambient_support(C::MixedCellCone)
    return C.ambientSupport
end

"""
    points(C::MixedCellCone)

Return the ambient points that define the support of the mixed cell cone `C`.
"""
function points(C::MixedCellCone)
    return points(ambient_support(C))
end

"""
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