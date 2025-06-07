###############################################################################
#
#  Chains of flats
#
###############################################################################

@doc raw"""
    ChainOfFlats

A chain of flats is an ascending sequence of flats of a matroid, starting from the empty set and ending at the ground set. Note that we do not include these two sets inside the vector in the `flats` field.
"""
struct ChainOfFlats
    matroid::Union{RealisableMatroid,Matroid}
    flats::Vector{Flat}
end



###############################################################################
#
#  Accessors
#
###############################################################################

@doc raw"""
    matroid(C::ChainOfFlats)

Return the matroid that `C` is a chain of flats of.
"""
function matroid(C::ChainOfFlats)
    return C.matroid
end

@doc raw"""
    flats(C::ChainOfFlats)

Return the vector of non trivial flats in `C`.
"""
function flats(C::ChainOfFlats)
    return C.flats
end

Base.getindex(C::ChainOfFlats, i::Int) = flats(C)[i]



###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    chain_of_flats(M::Union{RealisableMatroid, Matroid}, flats::Vector{Flat})

Construct a chain of flats from a matroid and a vector of flats.
"""
function chain_of_flats(M::Union{RealisableMatroid,Matroid}, flats::Vector{Flat})
    # check that we have a valid chain of flats
    # @assert !isempty(first(flats)) "First flat cannot be the empty set"
    # @assert !isequal(elements(last(flats)), ground_set(M)) "Last flat cannot be the ground set"
    # @assert is_subsequence([Set(elements(f)) for f in flats], Set.(Oscar.flats(matroid(M)))) "Did not provide a valid chain of flats"
    # @assert all(length(elements(flats[i])) < length(elements(flats[i+1])) for i in 1:length(flats)-1) "Flats must be strictly increasing in length"
    return ChainOfFlats(M, flats)
end

function chain_of_flats(M::Union{RealisableMatroid,Matroid}, flats::Vector{Vector{Int}})
    return chain_of_flats(M, [Flat(M, Set(f)) for f in flats])
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



###############################################################################
#
#  Conversions
#
###############################################################################

function Base.convert(::Set{Int}, F::Flat)
    return elements(F)
end



###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, C::ChainOfFlats)
    if isempty(flats(C))
        print(io, "∅ ⊊ {" * join(sort(collect(ground_set(matroid(C)))), ", ") * "}")
        return
    end
    # Convert each flat to a set and join with ⊊
    flat_strings = ["{" * join(sort(collect(elements(f))), ", ") * "}" for f in flats(C)]
    print(io, "∅ ⊊ " * join(flat_strings, " ⊊ ") * " ⊊ {" * join(sort(collect(ground_set(matroid(C)))), ", ") * "}")
end



###############################################################################
#
#  Properties
#
###############################################################################

@doc raw"""
    colength(C::ChainOfFlats)

Return the colength of the chain of flats `C`, i.e., the rank of its matroid minus its length.
"""
function colength(C::ChainOfFlats)
    return length(C) - rank(matroid(C)) + 1
end

@doc raw"""
    length(C::ChainOfFlats)

Return the length of the chain of flats `C`.
"""
function Base.length(C::ChainOfFlats)
    return length(flats(C))
end

@doc raw"""
    is_maximal(C::ChainOfFlats)

Return `true` if the chain of flats `C` is maximal and `false` otherwise.
"""
function is_maximal(C::ChainOfFlats)
    return iszero(colength(C))
end

function Base.isempty(C::ChainOfFlats)
    return isempty(flats(C))
end

function Base.isequal(C::ChainOfFlats, D::ChainOfFlats)
    return flats(C) == flats(D) && matroid(C) == matroid(D)
end

function Base.:(==)(C::ChainOfFlats, D::ChainOfFlats)
    return flats(C) == flats(D) && matroid(C) == matroid(D)
end

@doc raw"""
    Base.:(<)(C::ChainOfFlats, D::ChainOfFlats)

Check if the chain of flats `C` is a proper subsequence of the chain of flats `D`.
"""
function Base.:<(C::ChainOfFlats, D::ChainOfFlats)
    return (length(flats(C)) < length(flats(D))) && is_subsequence(flats(C), flats(D))
end

function is_subsequence(sub::Vector{T}, vec::Vector{T})::Bool where T
    # If sub is empty, it's technically a subsequence
    isempty(sub) && return true

    # Keep track of the last matched index in vec
    last_matched_index = 0

    # Iterate through each element in sub
    for s in sub
        # Find the index of s in vec, starting after the last matched index
        found_index = findnext(x -> is_equal(x, s), vec, last_matched_index + 1)

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
            if isnothing(findfirst(is_equal(v, entries(x)) for x in looplessFaceVertices))
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

    # add all ones vector


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
            push!(newFlats, elements(flats(C)[i]))
        else
            push!(newFlats, setdiff(elements(flats(C)[i]), elements(flats(C)[i-1])))
        end
    end

    push!(newFlats, setdiff(ground_set(matroid(C)), elements(flats(C)[end])))

    return newFlats
end


@doc raw"""
    maximal_refinements(C::ChainOfFlats)

Returns all maximal chains of flats that are refinements of `C`.

A refinement is a chain that contains the original chain as a subsequence.
Maximal chains are those that cannot be refined any further.
"""
function maximal_refinements(C::ChainOfFlats)::Vector{ChainOfFlats}
    mat = matroid(C)

    # Augment the chain with the empty set and the ground set.
    full_chain = [empty_flat(mat); flats(C); ground_flat(mat)]

    # Helper function: given two flats F and G (with F ⊂ G), return all flats F' with F ⊂ F' ⊂ G.
    function intermediate_flats(mat, F::Flat, G::Flat)::Vector{Flat}
        candidates = Set{Set{Int}}()
        for e in setdiff(elements(G), elements(F))
            # Compute the closure of F augmented by e.
            candidate = closure(mat, union(elements(F), [e]))
            # We want candidates strictly between F and G, and the candidate must be a flat.
            # (Since candidate = closure(candidate) by construction, it is a flat.)
            if elements(F) ⊊ candidate && candidate ⊊ elements(G)
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

function breaking_direction(maximalChainOfFlats::ChainOfFlats, nonmaximalChainOfFlats::ChainOfFlats)

    @assert matroid(maximalChainOfFlats) == matroid(nonmaximalChainOfFlats) "The matroids of the chains of flats must be the same"

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



# Given a maximal chain of flats and a subchain thereof of colength 1,
# return the star of the corresponding maximal Bergman cone around the facet of codimension 1
function bergman_cone_star(maximalChainOfFlats::ChainOfFlats, nonmaximalChainOfFlats::ChainOfFlats)
    println("maximalChainOfFlats: ", maximalChainOfFlats)
    println("nonmaximalChainOfFlats: ", nonmaximalChainOfFlats)
    @assert matroid(maximalChainOfFlats) == matroid(nonmaximalChainOfFlats) "The matroids of the chains of flats must be the same"
    @assert length(flats(maximalChainOfFlats)) == length(flats(nonmaximalChainOfFlats))+1 "nonmaximalChainOfFlats must have colength 1"

    facetRays = [ indicator_vector(Fj) for Fj in full_flats(nonmaximalChainOfFlats) if !isempty(Fj) ]

    println("flats(maximalChainOfFlats): ", flats(maximalChainOfFlats))
    println("flats(nonmaximalChainOfFlats): ", flats(nonmaximalChainOfFlats))

    extraFlatIndex = findfirst(Fj->!(Fj in flats(nonmaximalChainOfFlats)), flats(maximalChainOfFlats)) # TODO: this can be sped up
    @assert !isnothing(extraFlatIndex) "There must be a flat that is different in the two chains of flats"
    extraFlat = flats(maximalChainOfFlats)[extraFlatIndex]
    extraRay = indicator_vector(extraFlat)

    println("facetRays: ", facetRays)
    println("extraRay: ", extraRay)


    return positive_hull([extraRay],facetRays)
end


# Given a chain of flats, a subchain thereof, and a tiebreaking weight vector,
# return true if the subchain together with the tiebreaking vector induces the maximal chain of flats.
function induces_refinement(chain::ChainOfFlats, subchain::ChainOfFlats, tiebreaker::Vector)

    # basic input sanity checks
    @assert matroid(chain) == matroid(subchain) "The matroids of the chains of flats must be the same"
    @assert length(flats(chain)) == length(flats(subchain))+1 "nonmaximalChainOfFlats must have colength 1"

    setsInChain = TropicalHomotopyContinuation.elements.(full_flats(chain))

    # iterate over the reduced flats of the subchain, partition them further using the tiebreaker
    # and verify that all resulting flats are contained in the chain
    Finduced = Set{Int}()
    for Fred in reduced_flats(subchain)

        # partition Fred using the tiebreaker in the form of a Dict
        FredPartition = Dict{eltype(tiebreaker), Vector{Int}}()
        for i in Fred
            if !haskey(FredPartition, tiebreaker[i])
                FredPartition[tiebreaker[i]] = Int[]
            end
            push!(FredPartition[tiebreaker[i]], i)
        end

        # iterate over keys of the Dict in descending order,
        # and add them to Finduced and check if it is in the chain
        for key in sort(collect(keys(FredPartition)), rev=true)
            Finduced = union(Finduced, FredPartition[key])
            if !(Finduced in setsInChain)
                # induced flat not in chain
                return false
            end
        end
    end

    # all induced flats were in the chain
    return true
end
