using TropicalHomotopyContinuation
using Oscar

# specify the Bergman fan in the form of a realisation matrix
# for realisation numbers of minimally rigid graphs,
# it consists of two copies of the vertex-edge matrix
linearMatrix = matrix(QQ, [ 1  1  0  0  0  0  0  0  0  0;
                           -1  0 -1 -1  0  0  0  0  0  0;
                            0 -1  1  0 -1  0  0  0  0  0;
                            0  0  0  1  1  0  0  0  0  0;
                            0  0  0  0  0  1  1  0  0  0;
                            0  0  0  0  0 -1  0 -1 -1  0;
                            0  0  0  0  0  0 -1  1  0 -1;
                            0  0  0  0  0  0  0  0  1  1])
linearMatrix = echelon_form(linearMatrix, trim=true)
M = matroid(linearMatrix) # column matroid


# specify the target hypersurfaces in the form of tropical polynomials
# for realisation numbers, these are xi*yi+0 and x1+0 to cut away lineality
T = tropical_semiring()
Txy,x,y = polynomial_ring(T, "x"=>1:5, "y"=>1:5)
F = vcat([xi*yi+0 for (xi, yi) in zip(x, y)],[x[1]+0])
targetSupport = mixed_support(F)


# construct the starting mixed support and the starting mixed cells
Oscar.randseed!(13371337)
startingSupport, startingCell = starting_data(targetSupport, M)


# construct the tracker from the starting mixed support to the target mixed support
# possible paths are :coefficient_wise (fast, non-deterministic) and :straight_line (slow, deterministic)
# T = tracker(startingSupport, targetSupport, [startingCell], path=:coefficient_wise)
T = tracker(startingSupport, targetSupport, [startingCell], path=:straight_line)


# Move tracker until reaching endgame
AbstractAlgebra.set_verbosity_level(:TropicalHomotopyContinuation, 1)
@time track!(T)


# Perform endgame
@time endgame(T)



#################################################################################
#
#  Comparison with naive stable intersection
#
#################################################################################

# @time begin
#     toIntersect = tropical_hypersurface.(F)
#     toIntersect = vcat(toIntersect, [tropical_linear_space(linearMatrix)])
#     reduce(stable_intersection,toIntersect)
# end
