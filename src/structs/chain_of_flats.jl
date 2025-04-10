export Flat, ChainOfFlats, chain_of_flats, flat, closure, ground_set, empty_flat, ground_flat, is_subsequence, indicator_vector, maximal_refinements, breaking_direction

struct Flat
    matroid::Union{RealisableMatroid,Matroid}
    basis::Set{Int}
end

@doc raw"""
    ChainOfFlats

A chain of flats is an ascending sequence of flats of a matroid, starting from the empty set and ending at the ground set. Note that we do not include these two sets inside the vector in the `flats` field.
"""
struct ChainOfFlats
    matroid::Union{RealisableMatroid,Matroid}
    flats::Vector{Flat}
end

import Oscar.flats
@doc raw"""
    flats(C::ChainOfFlats)

Return the vector of non trivial flats in `C`.
"""
function flats(C::ChainOfFlats)
    return C.flats
end

@doc raw"""
    matroid(C::ChainOfFlats)

Return the matroid that `C` is a chain of flats of.
"""
function matroid(C::ChainOfFlats)
    return C.matroid
end

@doc raw"""
    matroid(F::Flat)

Return the matroid that `F` is a flat of.
"""
function matroid(F::Flat)
    return F.matroid
end

@doc raw"""
    basis(F::Flat)

Return the basis of a flat.
"""
function basis(F::Flat)
    return F.basis
end

@doc raw"""
    chain_of_flats(M::Union{RealisableMatroid, Matroid}, flats::Vector{Flat})

Construct a chain of flats from a matroid and a vector of flats.
"""
function chain_of_flats(M::Union{RealisableMatroid,Matroid}, flats::Vector{Flat})
    # check that we have a valid chain of flats
    if isempty(flats)
        return ChainOfFlats(M, Flat[])
    end
    # disable all checks for now
    # @assert !isempty(first(flats)) "First flat cannot be the empty set"
    # @assert !isequal(basis(last(flats)), ground_set(M)) "Last flat cannot be the ground set"
    # @assert is_subsequence([Set(basis(f)) for f in flats], Set.(Oscar.flats(matroid(M)))) "Did not provide a valid chain of flats"
    # @assert all(length(basis(flats[i])) < length(basis(flats[i+1])) for i in 1:length(flats)-1) "Flats must be strictly increasing in length"

    return ChainOfFlats(M, flats)
end

function chain_of_flats(M::Union{RealisableMatroid,Matroid}, flats::Vector{Vector{Int}})
    return chain_of_flats(M, [Flat(M, Set(f)) for f in flats])
end

@doc raw"""
    flat(M::Union{RealisableMatroid, Matroid}, basis::Union{Set{Int}, Set{}})

Construct a flat from a matroid and a basis. Raises an error if the basis is not a valid flat.
"""
function flat(M::Matroid, basis::Set{Int})
    # check that the elements of basis actually index a basis in M
    @assert basis in Set.(Oscar.flats(M)) "Did not provide a valid basis for a flat"
    # @assert !isequal(basis, ground_set(M)) "Basis cannot be the ground set"
    return Flat(M, basis)
end

@doc raw"""
    flat(M::Union{RealisableMatroid, Matroid}, basis::Union{Set{Int}, Set{}})

Construct a flat from a matroid and a basis. Raises an error if the basis is not a valid flat.
"""
function flat(M::RealisableMatroid, basis::Set{Int})
    @assert closure(M, basis) == basis "The basis is not a valid flat"
    return flat(matroid(M), basis)
end

@doc raw"""
    flat(M::Union{RealisableMatroid, Matroid}, basis::Union{Vector{Int}, Vector{}})

Construct a flat from a matroid and a basis. Raises an error if the basis is not a valid flat.
"""
function flat(M::Union{RealisableMatroid,Matroid}, basis::Vector{Int})
    return flat(M, Set(basis))
end

function Base.show(io::IO, F::Flat)
    print(io, "Flat with basis $(basis(F))")
end

function Base.show(io::IO, C::ChainOfFlats)
    if isempty(flats(C))
        print(io, "∅ ⊊ {" * join(sort(collect(ground_set(matroid(C)))), ", ") * "}")
        return
    end
    # Convert each basis to a set and join with ⊊
    flat_strings = ["{" * join(sort(collect(basis(f))), ", ") * "}" for f in flats(C)]
    print(io, "∅ ⊊ " * join(flat_strings, " ⊊ ") * " ⊊ {" * join(sort(collect(ground_set(matroid(C)))), ", ") * "}")
end

@doc raw"""
    is_subsequence(sub::Vector{T}, vec::Vector{T})::Bool where T

Check if `sub` is a subsequence of `vec`.
"""
function is_subsequence(sub::Vector{T}, vec::Vector{T})::Bool where T
    # If sub is empty, it's technically a subsequence
    isempty(sub) && return true

    # Keep track of the last matched index in vec
    last_matched_index = 0

    # Iterate through each element in sub
    for s in sub
        # Find the index of s in vec, starting after the last matched index
        found_index = findnext(x -> isequal(x, s), vec, last_matched_index + 1)

        # If not found, or found at an earlier index, return false
        if isnothing(found_index)
            return false
        end

        # Update the last matched index
        last_matched_index = found_index
    end

    return true
end

@doc raw"""
    ground_set(M::Union{RealisableMatroid, Matroid})

Return the ground set of a matroid as a set.
"""
function ground_set(M::Matroid)
    return Set(M.groundset)
end

@doc raw"""
    Base.length(C::ChainOfFlats)

Return the length of the chain of flats `C`.
"""
function Base.length(C::ChainOfFlats)
    return length(flats(C))
end

@doc raw"""
    is_maximal(C::ChainOfFlats)

Check if the chain of flats `C` is maximal.
"""
function is_maximal(C::ChainOfFlats)
    return length(C) == rank(matroid(C)) - 1
end

function Base.convert(::Set{Int}, F::Flat)
    return basis(F)
end

function Base.isempty(F::Flat)
    return isempty(basis(F))
end

function Base.isempty(C::ChainOfFlats)
    return isempty(flats(C))
end

function Base.isequal(F::Flat, G::Flat)
    return basis(F) == basis(G) && matroid(F) == matroid(G)
end

function Base.isequal(C::ChainOfFlats, D::ChainOfFlats)
    return flats(C) == flats(D) && matroid(C) == matroid(D)
end

@doc raw"""
    indicator_vector(flat::Flat)

Return the indicator vector of a flat.

This is a vector of length equal to the cardinality of the ground set of the underlying matroid of `flat`, with 1s at the indices in the basis of the flat, and 0s otherwise.
"""
function indicator_vector(flat::Flat)
    v = zeros(Int, length(ground_set(matroid(flat))))
    for i in basis(flat)
        v[i] = 1
    end
    return v
end

@doc raw"""
    Base.:(<)(C::ChainOfFlats, D::ChainOfFlats)

Check if the chain of flats `C` is a proper subsequence of the chain of flats `D`.
"""
function Base.:<(C::ChainOfFlats, D::ChainOfFlats)
    return (length(flats(C)) < length(flats(D))) && is_subsequence(flats(C), flats(D))
end

@doc raw"""
    chain_of_flats(M::Union{RealisableMatroid, Matroid}, w::TropicalPoint)

Construct a chain of flats induced on the matroid `M` by the point `w` in the Bergman fan of `M`.

Note that is required that the length of `w` is equal to the size of the ground set of `M`.
"""
function chain_of_flats(M::Union{RealisableMatroid,Matroid}, w::TropicalPoint)
    @assert length(w) == length(ground_set(M)) "The tropical point must have the same length as the ground set of the matroid"

    w = -w
    # Create a dictionary to group indices by their values
    value_indices = Dict{eltype(w),Vector{Int}}()

    # Populate the dictionary
    for (idx, val) in enumerate(w)
        if !haskey(value_indices, val)
            value_indices[val] = Int[]
        end
        push!(value_indices[val], idx)
    end

    # Sort the grouped indices
    sorted_groups = [sort(indices) for (_, indices) in sort(collect(value_indices), by=x -> x[1])]

    # Create cumulative union
    flat_indices = Vector{Int}[]
    cumulative_union = Int[]

    for group in sorted_groups[1:end-1]  # Exclude the last group
        cumulative_union = sort(union(cumulative_union, group))
        push!(flat_indices, cumulative_union)
    end

    return chain_of_flats(M, flat_indices)

end

function full_flats(C::ChainOfFlats)
    return [empty_flat(matroid(C)); flats(C); ground_flat(matroid(C))]
end

@doc raw"""
    loopless_face(C::ChainOfFlats)

Return the vertices of the loopless face whose Bergman cone contains the cone dual to `C`.
"""
function loopless_face(C::ChainOfFlats)
    looplessFaceVertices = Point[]

    candidateBases = Set.(vec(collect(Iterators.product(reduced_flats(C)...))))


    for candidateBasis in candidateBases
        if is_basis(matroid(C), candidateBasis)
            # Create indicator vector for the candidate basis
            v = zeros(Int, length(ground_set(matroid(C))))
            for i in candidateBasis
                v[i] = -1
            end
            if isnothing(findfirst(isequal(v, entries(x)) for x in looplessFaceVertices))
                push!(looplessFaceVertices, Point(v))
            end
        end
    end

    return looplessFaceVertices
end

@doc raw"""
    cone(C::ChainOfFlats)

Return the fine structure cone dual to the chain of flats `C`.
"""
function cone(C::ChainOfFlats)

    reducedFlats = reduced_flats(C)

    equalities = Vector{QQFieldElem}[]
    inequalities = Vector{QQFieldElem}[]

    for (i, F) in enumerate(reducedFlats)
        F1, Frest = Iterators.peel(F)
        for Fj in Frest
            equality = zeros(QQ, length(ground_set(matroid(C))))
            equality[F1] = 1
            equality[Fj] = -1
            push!(equalities, equality)
        end
        
        for j in 1:(i-1)
            G = reducedFlats[j]
            for g in G
                inequality = zeros(QQ, length(ground_set(matroid(C))))
                inequality[g] = -1
                inequality[F1] = 1
                push!(inequalities, inequality)
            end
        end
        
    end

    @debug "Equalities: $(equalities)"
    @debug "Inequalities: $(inequalities)"

    # if length(equalities) == 0
    #     return cone_from_inequalities(Oscar.matrix(QQ, inequalities))
    # else
    #     return cone_from_inequalities(Oscar.matrix(QQ, inequalities), Oscar.matrix(QQ, equalities))
    # end

    return inequalities, equalities

end

@doc raw"""
    reduced_flats(C::ChainOfFlats)

Return the reduced flats of the chain of flats `C`. This is a list of sets of indices, where each set is the indices of the elements in the flat that are not in the previous flat.
"""
function reduced_flats(C::ChainOfFlats)

    if isempty(flats(C))
        return Set{Int}[]
    end
    
    newFlats = Set{Int}[]

    for i in 1:length(flats(C))
        if i == 1
            push!(newFlats, basis(flats(C)[i]))
        else
            push!(newFlats, setdiff(basis(flats(C)[i]), basis(flats(C)[i-1])))
        end
    end

    push!(newFlats, setdiff(ground_set(matroid(C)), basis(flats(C)[end])))

    return newFlats
end

@doc raw"""
    ChainOfFlatsCone

A chain of flats cone is a cone defined by a chain of flats, along with a set of equations and inequalities that define the fine structure cone.
"""
struct ChainOfFlatsCone

    chainOfFlats::ChainOfFlats
    equations::Vector{Vector{QQFieldElem}}
    inequalities::Vector{Vector{QQFieldElem}}

end

function chain_of_flats_cone(chainOfFlats::ChainOfFlats, equations::Vector{Vector{QQFieldElem}}, inequalities::Vector{Vector{QQFieldElem}})

    # Check that the equations and inequalities are of the same length as the ground set of the matroid
    @assert all(length(equations[i]) == length(ground_set(matroid(chainOfFlats))) for i in 1:length(equations)) "The equations and inequalities must be of the same length as the ground set of the matroid"

    return ChainOfFlatsCone(chainOfFlats, equations, inequalities)

end

function polymake_cone(C::ChainOfFlatsCone)::Cone

    A = Oscar.matrix(QQ, C.inequalities)
    b = Oscar.matrix(QQ, C.equations)

    return cone_from_inequalities(A, b)

end

function chain_of_flats(C::ChainOfFlatsCone)
    return C.chainOfFlats
end

function inequalities(C::ChainOfFlatsCone)
    return C.inequalities
end

function Base.show(io::IO, C::ChainOfFlatsCone)

    print(io, "Cone defined by the chain of flats $(chain_of_flats(C))")

end

function Base.in(w::TropicalPoint, C::ChainOfFlatsCone)

    return all([dot(v, w) <= 0 for v in inequalities(C)]) && all([dot(v, w) == 0 for v in equations(C)])

end

function empty_flat(M::Union{RealisableMatroid,Matroid})
    return Flat(M, Set{Int}())
end

function ground_flat(M::Union{RealisableMatroid,Matroid})
    return Flat(M, ground_set(M))
end

function Base.:(==)(F::Flat, G::Flat)
    return basis(F) == basis(G) && matroid(F) == matroid(G)
end

function Base.:(==)(C::ChainOfFlats, D::ChainOfFlats)
    return flats(C) == flats(D) && matroid(C) == matroid(D)
end

Base.getindex(C::ChainOfFlats, i::Int) = flats(C)[i]

"""
    maximal_refinements(C::ChainOfFlats)

Returns all maximal chains of flats that are refinements of `C`.

A refinement is a chain that contains the original chain as a subsequence.
Maximal chains are those that cannot be refined any further.
"""
function maximal_refinements(C::ChainOfFlats)::Vector{ChainOfFlats}
    mat = matroid(C)

    # Augment the chain with the empty set and the ground set.
    # (Assuming you have functions `empty_flat(mat)` and `ground_flat(mat)`.)
    full_chain = [empty_flat(mat); flats(C); ground_flat(mat)]

    # Helper function: given two flats F and G (with F ⊂ G), return all flats F' with F ⊂ F' ⊂ G.
    function intermediate_flats(mat, F::Flat, G::Flat)::Vector{Flat}
        candidates = Set{Set{Int}}()
        # Consider each element in G.basis that is not already in F.basis.
        for e in setdiff(basis(G), basis(F))
            # Compute the closure of F augmented by e.
            candidate = closure(mat, union(basis(F), [e]))
            # We want candidates strictly between F and G, and the candidate must be a flat.
            # (Since candidate = closure(candidate) by construction, it is a flat.)
            if basis(F) ⊊ candidate && candidate ⊊ basis(G)
                push!(candidates, Set(candidate))
            end
        end
        return [Flat(mat, s) for s in candidates]
    end

    # Recursive helper: refine the chain starting at gap i (i.e. between full_chain[i] and full_chain[i+1]).
    function refine_chain(chain::Vector{Flat}, i::Int)::Vector{Vector{Flat}}
        # If we've reached the end of the chain, return the chain as is.
        if i == length(chain)
            return [chain]
        end

        # Find intermediate flats between chain[i] and chain[i+1]
        intermediates = intermediate_flats(mat, chain[i], chain[i+1])

        # If no refinement is possible in this gap, move to the next gap.
        if isempty(intermediates)
            return refine_chain(chain, i + 1)
        else
            refined_chains = Vector{Vector{Flat}}()
            # For each possible intermediate flat, insert it and then try to refine further.
            for f in intermediates
                # Create a new chain by inserting f at position i+1
                new_chain = copy(chain)
                insert!(new_chain, i+1, f)
                # Recursively refine from the same gap (as more flats might be inserted in the new gap)
                for refined in refine_chain(new_chain, i+1)
                    push!(refined_chains, refined)
                end
            end
            return refined_chains
        end
    end

    # Start the recursive refinement from the first gap.
    all_full_chains = refine_chain(full_chain, 1)
    # Remove the initial empty and ground flats before returning.
    result = [chain_of_flats(mat, ch[2:end-1]) for ch in all_full_chains]
    return result
end

@doc raw"""
    closure(M::Union{RealisableMatroid, Matroid}, elements::Set{Int})

Compute the closure of a set of elements in a matroid.
"""
function closure(M::Matroid, elements::Set{Int})
    # Fallback implementation using flats
    # Find the smallest flat containing elements
    for f in Oscar.flats(M)
        if all(e in f for e in elements)
            return f
        end
    end
    return ground_set(M)  # Default if no flat is found
    
end

@doc raw"""
    closure(M::RealisableMatroid, elements::Set{Int})

Compute the closure of a set of elements in a realisable matroid.
"""
function closure(M::RealisableMatroid, elements::Set{Int})

    mat = matrix(M)

    # check that including any other element keeps the rank the same
    for i in 1:length(ground_set(M))
        if !(i in elements)
            if Oscar.rank(mat[:, [collect(elements); i]]) == Oscar.rank(mat[:, collect(elements)])
                push!(elements, i)
            end
        end
    end

    return elements
end

function breaking_direction(maximalChainOfFlats::ChainOfFlats, nonmaximalChainOfFlats::ChainOfFlats)

    M = matroid(maximalChainOfFlats)
    @assert M == matroid(nonmaximalChainOfFlats) "The matroids of the chains of flats must be the same"

    maximalChain = full_flats(maximalChainOfFlats)
    nonMaximalChain = full_flats(nonmaximalChainOfFlats)

    # find the first flat that is different
    changingFlat = 2
    for i in 2:length(maximalChain)-1
        if maximalChain[i] != nonMaximalChain[i]
            changingFlat = i
            break
        end
    end

    return 2*indicator_vector(maximalChain[changingFlat]) - indicator_vector(maximalChain[changingFlat+1]) - indicator_vector(maximalChain[changingFlat - 1])
end
