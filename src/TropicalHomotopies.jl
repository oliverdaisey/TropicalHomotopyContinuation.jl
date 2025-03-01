module TropicalHomotopies

using Oscar
include("structs/point.jl")
include("structs/weight.jl")
include("structs/support.jl")
include("structs/mixed_support.jl")
include("structs/cayley_embedding.jl")
include("structs/mixed_cell_cone.jl")

# Write your package code here.
pts = [point(0,0), point(1,0), point(0,1), point(1,1)]
wts = [weight(0), weight(1), weight(0), weight(2)]

mySupport = support(pts, wts)

pts = [point(0,0), point(2,0), point(2,2)]
wts = [weight(0), weight(0), weight(0)]

mySecondSupport = support(pts, wts)

display(mySupport)
display(mySecondSupport)

mixedSupport = mixed_support((mySupport, mySecondSupport))

display(mixedSupport)

cayley = cayley_embedding(mixedSupport)

println("Constructed cayley embedding")
display(matrix(cayley))

# submatrix = cayley[mySecondSupport]

end