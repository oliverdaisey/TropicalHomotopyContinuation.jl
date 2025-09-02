###############################################################################
#
#  Chains of flats
#
###############################################################################

@doc raw"""
    Flat

A closed subset of the ground set of a matroid.  Here, closed means adding any ground set element increases its rank.
"""
struct Flat
    matroid::Union{RealisableMatroid,Matroid}
    elements::Set{Int}
end



###############################################################################
#
#  Accessors
#
###############################################################################

@doc raw"""
    matroid(F::Flat)

Return the matroid that `F` is a flat of.
"""
function matroid(F::Flat)
    return F.matroid
end

@doc raw"""
    elements(F::Flat)

Return the elements of `F`.
"""
function elements(F::Flat)
    return F.elements
end



###############################################################################
#
#  Constructors
#
###############################################################################

@doc raw"""
    flat(M::Union{RealisableMatroid, Matroid}, elements::Union{Set{Int}, Set{}})

Construct a flat from a matroid and a set of elements. Raises an error if the elements do not define a valid flat.
"""
function flat(M::Matroid, elements::Set{Int})
    # check that the elements actually index a flat in M
    @assert elements in Set.(Oscar.flats(M)) "Did not provide valid elements for a flat"
    # @assert !is_equal(elements, ground_set(M)) "Cannot be the ground set"
    return Flat(M, elements)
end

@doc raw"""
    flat(M::Union{RealisableMatroid, Matroid}, elements::Union{Set{Int}, Set{}})

Construct a flat from a matroid and elements. Raises an error if the elements do not index a valid flat.
"""
function flat(M::RealisableMatroid, elements::Set{Int})
    @assert closure(M, elements) == elements "The elements $(elements) do not index a valid flat"
    return Flat(M, elements)
end

@doc raw"""
    flat(M::Union{RealisableMatroid, Matroid}, elements::Union{Vector{Int}, Vector{}})

Construct a flat from a matroid and a set of elements. Raises an error if the elements do not index a valid flat.
"""
function flat(M::Union{RealisableMatroid,Matroid}, elements::Vector{Int})
    return flat(M, Set(elements))
end

function empty_flat(M::Union{RealisableMatroid,Matroid})
    return Flat(M, Set{Int}())
end

function ground_flat(M::Union{RealisableMatroid,Matroid})
    return Flat(M, ground_set(M))
end

@doc raw"""
    closure(M::RealisableMatroid, elements::Set{Int})

Compute the closure of a set of elements in a realisable matroid.
"""
function closure(M::RealisableMatroid, elems::Set{Int})

    # check if including any other element keeps the rank the same
    for i in setdiff(ground_set(M), elems)
        if rank(M, union(elems, Set{Int}([i]))) == rank(M, elems)
            push!(elems, i)
        end
    end

    return elems
end



###############################################################################
#
#  Printing
#
###############################################################################

function Base.show(io::IO, F::Flat)
    print(io, "Flat indexed by $(elements(F))")
end



###############################################################################
#
#  Properties
#
###############################################################################

function Base.isempty(F::Flat)
    return isempty(elements(F))
end

function isequal(F::Flat, G::Flat)
    return elements(F) == elements(G) && matroid(F) == matroid(G)
end

function Base.:(==)(F::Flat, G::Flat)
    return isequal(F,G)
end

function Base.hash(F::Flat, h::UInt)
    # Combine hashes of sorted elements to make sure set equality aligns with hash equality
    elements_hash = hash(sort(collect(F.elements)), h)
    matroid_hash = hash(F.matroid, h)
    return hash((matroid_hash, elements_hash), h)
end



###############################################################################
#
#  Other functions
#
###############################################################################

@doc raw"""
    indicator_vector(flat::Flat)

Return the indicator vector of a flat.

This is a vector of length equal to the cardinality of the ground set of the underlying matroid of `flat`, with 1s at the indices in the elements of the flat, and 0s otherwise.
"""
function indicator_vector(flat::Flat)
    v = zeros(Int, length(ground_set(matroid(flat))))
    for i in elements(flat)
        v[i] = 1
    end
    return v
end

@doc raw"""
    rank(F::Flat)

Return the rank of the flat `F`.
"""
function rank(F::Flat)
    return rank(matroid(F), elements(F))
end
