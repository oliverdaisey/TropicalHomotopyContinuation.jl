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
