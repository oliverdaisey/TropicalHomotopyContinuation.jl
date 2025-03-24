M = matroid(Oscar.matrix(QQ, [-1 -2 -3 -4; 11 13 15 19]))
TT = tropical_semiring()
R, (x,y,z,w) = polynomial_ring(TT, ["x","y","z","w"])

f = 1*x + y + (-1)*z + (-2)*w
g = 11*x + 13*y + 15*z + 19*w

A, b = find_tropical_point((f,g), M)

using Oscar
solve(A, b)