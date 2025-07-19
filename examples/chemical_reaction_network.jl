using TropicalHomotopyContinuation
using Oscar

# specify the Bergman fan in the form of a realisation matrix
# for realisation numbers of minimally rigid graphs,
# it consists of two copies of the vertex-edge matrix
linearEquationMatrix = matrix(QQ, [0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0   0   0   0   0   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   1   1   0   0   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   1   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1   1   1   1   0   0   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   1   1   0   0   0   0   0;
                           0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   1   1   0   0   0;
                                   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   1]);
linearRealisationMatrix = transpose(kernel(linearEquationMatrix, side=:right))
linearMatrix = echelon_form(linearRealisationMatrix, trim=true)
M = matroid(linearRealisationMatrix) # column matroid


# specify the target hypersurfaces in the form of tropical polynomials
# for realisation numbers, these are xi*yi+0 and x1+0 to cut away lineality
T = tropical_semiring()
Txy,x,y = polynomial_ring(T, "x"=>1:19, "y"=>1:53)
F = [x[1] + y[1],
     x[2] + y[2],
     x[2]*x[4] + y[3],
     x[2] + y[4],
     x[3] + y[5],
     x[14] + y[6],
     x[14] + y[7],
     x[5]*x[8] + y[8],
     x[5] + y[9],
     x[7] + y[10],
     x[16] + y[11],
     x[3]*x[6] + y[12],
     x[6]*x[11] + y[13],
     x[15] + y[14],
     x[17] + y[15],
     x[19] + y[16],
     x[19] + y[17],
     x[7]*x[9] + y[18],
     x[15] + y[19],
     x[17] + y[20],
     x[4]*x[10] + y[21],
     x[10] + y[22],
     x[10] + y[23],
     x[11] + y[24],
     x[18] + y[25],
     y[26] + 0,
     x[11]*x[12] + y[27],
     x[11] + y[28],
     x[13] + y[29],
     x[16] + y[30],
     x[18] + y[31],
     y[32] + 0,
     x[1] + y[33],
     x[2] + y[34],
     x[3] + y[35],
     x[14] + y[36],
     x[15] + y[37],
     y[38] + 0,
     x[4] + y[39],
     x[5] + y[40],
     x[6] + y[41],
     x[7] + y[42],
     x[16] + y[43],
     x[17] + y[44],
     x[18] + y[45],
     x[19] + y[46],
     y[47] + 0,
     x[8] + y[48],
     y[49] + 0,
     x[9] + y[50],
     y[51] + 0,
     x[12] + y[52],
     x[13] + y[53]]
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
