using TropicalHomotopies
using Oscar

# define matrix encoding linear ideal
linearMatrix = matrix(QQ, [1 1 0 0 0 0 0 0 0 0; -1 0 -1 -1 0 0 0 0 0 0; 0 -1 1 0 -1 0 0 0 0 0; 0 0 0 1 1 0 0 0 0 0; 0 0 0 0 0 1 1 0 0 0; 0 0 0 0 0 -1 0 -1 -1 0; 0 0 0 0 0 0 -1 1 0 -1; 0 0 0 0 0 0 0 0 1 1])

# reduce linearMatrix to row echelon form
linearMatrix = echelon_form(linearMatrix, trim=true)
R, (x1,x2,x3,x4,x5,x6,x7,x8,x9,x10) = polynomial_ring(QQ, 10)
M = matroid(linearMatrix)

# define hypersurface supports
# Oscar.randseed!(31415296) # seed to reproduce bug (2 mixed cells)
# Oscar.randseed!(127) # 3 mixed cells at the end
Oscar.randseed!(143)
targetSupports = Support[]
for i in [1,2,3,4,5]
        pi = TropicalHomotopies.point([n in [i,i+5] ? 1 : 0 for n in 1:10])
        p0 = TropicalHomotopies.point([0 for n in 1:10])
            push!(targetSupports, TropicalHomotopies.support([pi,p0],[0,0]))
end
# add an extra hypersurface to cut the common lineality space
push!(targetSupports, TropicalHomotopies.support([TropicalHomotopies.point([1,0,0,0,0,0,0,0,0,0]),TropicalHomotopies.point([0,0,0,0,0,0,0,0,0,0])],[1,3]))

# define the target support
targetSupport = mixed_support(targetSupports)
Δ, σ = starting_data(targetSupport, M)
# T = tracker(Δ, targetSupport, [σ], path=:coefficient_wise)
T = tracker(Δ, targetSupport, [σ], path=:straight_line)

# compute the stable intersection
AbstractAlgebra.set_verbosity_level(:TropicalHomotopies, 0)
@time track!(T)
