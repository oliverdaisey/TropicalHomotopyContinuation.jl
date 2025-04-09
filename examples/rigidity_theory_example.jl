using TropicalHomotopies
import Oscar.matrix, Oscar.echelon_form, Oscar.QQ, Oscar.polynomial_ring

targetSupports = Support[]
linearMatrix = matrix(QQ, [1 1 0 0 0 0 0 0 0 0; -1 0 -1 -1 0 0 0 0 0 0; 0 -1 1 0 -1 0 0 0 0 0; 0 0 0 1 1 0 0 0 0 0; 0 0 0 0 0 1 1 0 0 0; 0 0 0 0 0 -1 0 -1 -1 0; 0 0 0 0 0 0 -1 1 0 -1; 0 0 0 0 0 0 0 0 1 1])

# reduce linearMatrix to row echelon form
linearMatrix = echelon_form(linearMatrix, trim=true)
R, (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = polynomial_ring(QQ, 10)

M = matroid(linearMatrix)
for i in [1,2,3,4,5]
        pi = point([n in [i,i+5] ? 1 : 0 for n in 1:10])
        p0 = point([0 for n in 1:10])
            push!(targetSupports, support([pi,p0],[0,0]))
end

push!(targetSupports, support([point([1,0,0,0,0,0,0,0,0,0]),point([0,0,0,0,0,0,0,0,0,0])],[1,3]))

targetSupport = mixed_support(targetSupports)

Δ, σ = starting_data(targetSupport, M)

T = tracker(Δ, targetSupport, [σ])

# compute the stable intersection
using Random
Random.seed!(1)
@profview display(stable_intersection(T))