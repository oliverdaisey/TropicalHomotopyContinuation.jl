using Oscar
using Revise

# alias for a tuple of tropical polynmomials
TropicalTuple{N} = NTuple{N, AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}}
TropicalPoly = AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}} # alias for a tropical polynomial
PluckerVector = Pair{Vector{Vector{Int}}, Vector{Int}}

T = tropical_semiring()

function compute_starting_data(F::TropicalTuple, variables::Vector{TropicalPoly})

    k = length(F)
    n = length(variables) - 1

    # step 1: compute p_start

    # compute plueker indices
    pluecker_indices = subsets(collect(1:(n+1)), k+1)
    pluecker_vector = [0 for i in 1:length(pluecker_indices)]
    p_start = PluckerVector(pluecker_indices, pluecker_vector)


    # step 2: compute linear forms
    linear_forms = Vector{TropicalPoly}()
    for i in 1:k
        
        l_i = T(0)
        # compute coefficients
        c::Vector{Oscar.TropicalSemiringElem{typeof(min)}} = [T.(0) for j in 0:n]
        for j in 0:n
            if j > k
                c[j+1] = T(-1)
            elseif j == i
                c[j+1] = T(1)
            end

            l_i += c[j+1] * variables[j+1]
        end

        push!(linear_forms, l_i)

    end

    # step 3: raise linear forms to appropriate powers
    deg = [get_degree(F[i]) for i in 1:k]
    for i in 1:k
        linear_forms[i] = linear_forms[i]^deg[i] # sometimes you will get a "characteristic not known" error
    end

    # step 4: compute mixed cell
    Σ = Vector{Polyhedron}()

    for i in 1:k
        # I need to get better at Julia vector/matrix ops
        ei = [0 for j in 0:n]
        e0 = [0 for j in 0:n]
        e0[1] = 1
        ei[i+1] = 1
        M = matrix(QQ, vcat(e0', ei'))
        σ_i = deg[i] * convex_hull(M)
        push!(Σ, σ_i)
    end

    # construct vertices of σ_p
    M = matrix(QQ, zeros(QQ, n-k+1, n+1))
    M[1,1] = 1
    for i in 1:n-k
        M[i+1, k+1+i] = 1
    end

    σ_p = convex_hull(M)

    S = sum(σ_i for σ_i in Σ) + σ_p

    return (p_start, linear_forms, S)

end

# f::MPoly{T} where T <: RingElement
function get_degree(f)::Int

    exponent_matrix = f.exps # the vectors are the columns of this matrix, I want to iterate over them
    degree = 0
    for i in 1:ncols(exponent_matrix)
        monomial_degree = sum(exponent_matrix[:,i])
        if monomial_degree > degree
            degree = monomial_degree
        end
    end

    return Int(degree)
end

function run_example()
    R, (w, x, y, z) = T["w", "x", "y", "z"]
    f0 = w^2 + w*x + x^2 + w*y + x*y + y^2 + w*z + x*z + y*z + z^2
    variables = [w, x, y, z]

    starting_data = compute_starting_data((f0,), variables)

    return starting_data
end