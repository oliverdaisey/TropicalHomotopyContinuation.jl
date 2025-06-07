# This tests the same example coming from rigidity theory as the one in rigidity_theory_example.jl.

using Oscar

A = matrix(QQ,[1 1 0 0 0; -1 0 -1 -1 0; 0 -1 1 0 -1; 0 0 0 1 1])
M = vcat(hcat(A,zero_matrix(QQ,4,5)),hcat(zero_matrix(QQ,4,5),A))
TropL = tropical_linear_space(M)

TT = tropical_semiring()
R,(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5) = TT["x1","x2","x3","x4","x5","y1","y2","y3","y4","y5"]
TropH1 = tropical_hypersurface(x1*y1+0)
TropH2 = tropical_hypersurface(x2*y2+0)
TropH3 = tropical_hypersurface(x3*y3+0)
TropH4 = tropical_hypersurface(x4*y4+0)
TropH5 = tropical_hypersurface(x5*y5+0)

@time TropV = reduce(Oscar.stable_intersection,[TropH1,TropH2,TropH3,TropH4,TropH5,TropL])