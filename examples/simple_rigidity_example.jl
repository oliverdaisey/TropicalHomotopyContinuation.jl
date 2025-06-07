using TropicalHomotopyContinuation
using Oscar

A = matrix(QQ,[1 1 0 0 0; -1 0 -1 -1 0; 0 -1 1 0 -1; 0 0 0 1 1])
M = vcat(hcat(A,zero_matrix(QQ,4,5)),hcat(zero_matrix(QQ,4,5),A))

Oscar.randseed!(143)

TT = tropical_semiring()
R,(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) = TT["x1","x2","x3","x4","x5","y1","y2","y3","y4","y5"]

f1 = x1*y1 + 0
f2 = x2*y2 + 0
f3 = x3*y3 + 0
f4 = x4*y4 + 0
f5 = x5*y5 + 0

f0 = 1*x1 + 3 # to cut out lineality space

@time stable_intersection([f0,f1,f2,f3,f4,f5], matroid(M))