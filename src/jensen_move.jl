using Oscar

function jensen_time(T::Tracker, σ::MixedCell)::Union{QQFieldElem,Nothing}
    
    hypersurfaceDuals = transform_linear_support(chain_of_flats(σ))

    δ = combine(active_support(σ), hypersurfaceDuals)
    Δ = combine(ambient_support(T), hypersurfaceDuals)

    C = mixed_cell_cone(δ, Δ)
    
    @assert δ in C "The mixed cell being tracked is not in the mixed cell cone."

    v = direction(T)

    return minimum([dot(v, u) != 0 ? -dot(v, δ) / dot(v, u) : PosInf for u in facets(C)])
    
end
