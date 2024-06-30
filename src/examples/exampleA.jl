include("../main.jl")

function total_degree_starting_data(p::PluckerVector, F::TropicalTuple)

    k = length(F)
    n = ngens(parent(F[1])) - 1
    variables = gens(parent(F[1]))

    # step 1: compute p_start

    # compute plueker indices
    pluecker_indices = subsets(collect(1:(n+1)), k+1)
    pluecker_vector = [T(0) for i in 1:length(pluecker_indices)]
    p_start = PluckerVector(pluecker_indices, pluecker_vector)


    # step 2: compute linear forms
    F_start = Vector{TropicalPoly}()
    for i in 1:k
        
        l_i = zero(T)
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

        push!(F_start, l_i)

    end

    # step 3: raise linear forms to appropriate powers
    deg = [get_degree(F[i]) for i in 1:k]
    for i in 1:k
        F_start[i] = F_start[i]^deg[i] # sometimes you will get a "characteristic not known" error
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
    
    return (p_start, F_start, S)

end