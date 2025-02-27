module TropicalHomotopies

using Oscar
include("structs/point.jl")
include("structs/weight.jl")
include("structs/support.jl")
include("structs/mixed_support.jl")

# Write your package code here.
pts = [point(0,0), point(1,0), point(0,1), point(1,1)]
wts = [weight(0), weight(1), weight(0), weight(2)]

mySupport = support(pts, wts)

pts = [point(0,0), point(2,0), point(2,2)]
wts = [weight(0), weight(0), weight(0)]

mySecondSupport = support(pts, wts)

mixedSupport = mixed_support((mySupport, mySecondSupport))

export mixedSupport

end