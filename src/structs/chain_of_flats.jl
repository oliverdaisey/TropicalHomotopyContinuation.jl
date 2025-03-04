using Oscar

struct Flat
    matroid::Matroid
    basis::Set{Int}
end

@doc raw"""
    ChainOfFlats

A chain of flats is an ascending sequence of flats of a matroid, starting from the empty set and ending at the ground set.
"""
struct ChainOfFlats
    matroid::Matroid
    flats::Vector{Flat}
end

import Oscar.flats
@doc raw"""
    flats(C::ChainOfFlats)

Return the vector of flats of `C`.
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
    basis(F::Flat)

Return the basis of a flat.
"""
function basis(F::Flat)
    return F.basis
end

@doc raw"""
    chain_of_flats(M::Matroid, flats::Vector{Flat})

Construct a chain of flats from a matroid and a vector of flats.
"""
function chain_of_flats(M::Matroid, flats::Vector{Flat})
    # check that we have a valid chain of flats
    @assert !isempty(first(flats)) "First flat cannot be the empty set"
    @assert !isequal(basis(last(flats)), ground_set(M)) "Last flat cannot be the ground set"
    @assert is_subsequence([Set(basis(f)) for f in flats], Set.(Oscar.flats(M))) "Did not provide a valid chain of flats"
    @assert all(length(basis(flats[i])) < length(basis(flats[i+1])) for i in 1:length(flats)-1) "Flats must be strictly increasing in length"
    
    return ChainOfFlats(M, flats)
end

function chain_of_flats(M::Matroid, flats::Vector{Vector{Int}})
    return chain_of_flats(M, [Flat(M, Set(f)) for f in flats])
end

@doc raw"""
    flat(M::Matroid, basis::Union{Set{Int}, Set{}})

Construct a flat from a matroid and a basis. Raises an error if the basis is not a valid flat.
"""
function flat(M::Matroid, basis::Set{Int})
    # check that the elements of basis actually index a basis in M
    @assert basis in Set.(flats(M)) "Did not provide a valid basis for a flat"
    @assert !isequal(basis, ground_set(M)) "Basis cannot be the ground set"
    return Flat(M, basis)
end

@doc raw"""
    flat(M::Matroid, basis::Union{Vector{Int}, Vector{}})

Construct a flat from a matroid and a basis. Raises an error if the basis is not a valid flat.
"""
function flat(M::Matroid, basis::Vector{Int})
    return flat(M, Set(basis))
end

function Base.show(io::IO, F::Flat)
    print(io, "Flat with basis $(basis(F))")
end

function Base.show(io::IO, C::ChainOfFlats)
    print(io, "Chain of flats of length $(length(flats(C)))")
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
        found_index = findnext(x -> x == s, vec, last_matched_index + 1)
        
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
    ground_set(M::Matroid)

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