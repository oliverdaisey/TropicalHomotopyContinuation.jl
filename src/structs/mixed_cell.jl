@doc raw"""
    MixedCellCandidate

A candidate for a mixed cell, consisting of an active support and a chain of flats.
"""
struct MixedCell

    activeSupport::MixedSupport
    chainOfFlats::ChainOfFlats

end

@doc raw"""
    mixed_cell(activeSupport::MixedSupport, chainOfFlats::ChainOfFlats)::MixedCellCandidate

Construct a mixed cell candidate.
"""
function mixed_cell(activeSupport::MixedSupport, chainOfFlats::ChainOfFlats)::MixedCell
    return MixedCell(activeSupport, chainOfFlats)
end

@doc raw"""
    points(σ::MixedCell)

Return the points of the active support of the mixed cell candidate `σ`.
"""
function points(σ::MixedCell)
    return points(active_support(σ))
end

@doc raw"""
    chain_of_flats(σ::MixedCell)

Return the chain of flats of the mixed cell candidate `σ`.
"""
function chain_of_flats(σ::MixedCell)
    return σ.chainOfFlats
end

@doc raw"""
    active_support(σ::MixedCell)

Return the active support of the mixed cell candidate `σ`.
"""
function active_support(σ::MixedCell)
    return σ.activeSupport
end

@doc raw"""
    supports(σ::MixedCell)

Return the supports defining the active support of the mixed cell candidate `σ`.
"""
function supports(σ::MixedCell)
    return supports(active_support(σ))
end

function Base.length(σ::MixedCell)
    return length(supports(σ))
end

@doc raw"""
    Base.getindex(σ::MixedCell, p::Point)

Do not index mixed cell candidates directly as their data may be stale. Index the Tracker instead.
"""
function Base.getindex(σ::MixedCell, p::Point)
    error("Do not index mixed cells directly as their data may be stale. Index the Tracker instead.")
end

function Base.show(io::IO, σ::MixedCell)
    print(io, "Mixed cell with active points $(join(["$(p) ∈ S_$(findfirst(x -> p in x, supports(σ)))" for p in sort(points(σ))], ", ")) and chain of flats $(chain_of_flats(σ))")
end

@doc raw"""
    transform_linear_support(σ::MixedCell)

Returns the mixed support obtained by replacing the chain of flats with a tuple of hypersurfaces whose intersection equals the span of the corresponding fine structure cone.
"""
function transform_linear_support(C::ChainOfFlats)::MixedSupport

    gens = indicator_vector.(flats(C))
    push!(gens, ones(Int, length(ground_set(matroid(C)))))
    # find basis of orthogonal complement of span(gens)
    gensMatrix = Oscar.matrix(gens)
    _, kernel = Oscar.nullspace(gensMatrix)
    # convert kernel to list of vectors
    kernel = [kernel[:, i] for i in 1:ncols(kernel)]
    # add the vectors to pts
    pts = [point(g) for g in kernel]
    S = Support[]
    for p in pts
        push!(S, support([point(zeros(Int, length(p))), p], [0, 0]))
    end


    return mixed_support(tuple(S...))
end

@doc raw"""
    swap(σ::MixedCell, p::Point, q::Point)

Swap the point `p` for the new point `q` in the active support of the mixed cell candidate `σ`.
"""
function swap(σ::MixedCell, p::Point, q::Point)
    Δ = active_support(σ)
    S = supports(σ)
    index = findfirst(x -> p in x, S)
    # remove p from this support
    newSupport = support(S[index], q)
    # remove the key corresponding to p
    delete!(newSupport.entries, p)
    newSupports = [S[i] for i in 1:length(S)]
    newSupports[index] = newSupport
    return mixed_cell(mixed_support(tuple(newSupports...)), chain_of_flats(σ))
end

function has_same_active_support(σ::MixedCell, τ::MixedCell)
    if length(supports(σ)) != length(supports(τ))
        return false
    end

    # each support should be a subset of another support
    for S in supports(σ)
        if !any([is_subset(S, T) for T in supports(τ)])
            return false
        end
    end
    for T in supports(τ)
        if !any([is_subset(T, S) for S in supports(σ)])
            return false
        end
    end

    return true

end

function is_transverse(σ::MixedCell)

    return true

    Δ = active_support(σ)

    rows = Vector{Int}[]

    # add in rows from chain of flats
    pts = loopless_face(chain_of_flats(σ))
    p1 = first(pts)

    for p in pts
        if !is_equal(p1, p)
            push!(rows, p1 - p)
        end
    end

    for S in supports(Δ)
        p1 = first(points(S))
        for p in points(S)
            if !is_equal(p1, p)
                push!(rows, p1 - p)
            end
        end
    end

    flag, _ = Oscar.is_invertible_with_inverse(Oscar.matrix(QQ, rows))

    return flag

end

function are_support_heights_finite(Δ::MixedSupport, σ::MixedCell)
    for p in points(active_support(σ))
        if isinf(Δ[p])
            return false
        end
    end

    return true

end

function tropical_polyhedra_spans_and_mults(σ::MixedCell)

    polyhedra = []
    mults = []
    for S in supports(active_support(σ))

        # compute multiplicity
        push!(mults, length(lattice_points(convex_hull(entries.(points(S)))))-1)


        rows = Vector{QQFieldElem}[]
        pts = points(S)
        p1 = first(pts)
        for p in pts
            if !is_equal(p1, p)
                push!(rows, convert(Vector{QQFieldElem}, p1 - p))
            end
        end
        M = Oscar.matrix(QQ, rows)
        A = transpose(kernel(M, side=:right))
        push!(polyhedra, convex_hull(zero_matrix(QQ, 1, ncols(A)), zero_matrix(QQ, 0, ncols(A)), A))

    end

    d = max(ambient_dim.(polyhedra)...)

    # add columns for each indicator vector of chain
    cols = Vector{QQFieldElem}[]
    push!(cols, indicator_vector.(full_flats(chain_of_flats(σ)))...)
    # remove all zero vector from cols
    cols = [col for col in cols if col != zeros(QQ, length(col))]
    A = Oscar.matrix(QQ, cols)


    push!(polyhedra, convex_hull(zero_matrix(QQ, 1, ncols(A)), zero_matrix(QQ, 0, ncols(A)), A))
    push!(mults, 1)

    return polyhedra, mults

end

# slow
function multiplicity(σ::MixedCell)

    polyhedra, mults = tropical_polyhedra_spans_and_mults(σ)

    (C, m), rest = Iterators.peel(zip(polyhedra, mults))
    for (D, n) in rest
        m *= n*Oscar.tropical_intersection_multiplicity(C, D)
        C = intersect(C, D)
    end

    return m

end
