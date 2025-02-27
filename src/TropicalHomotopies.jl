module TropicalHomotopies

using Oscar
include("structs/point.jl")
include("structs/weight.jl")
include("structs/support.jl")
include("structs/mixed_subdivision.jl")

# Write your package code here.
myPoint = point(3,4)
mySecondPoint = point(5,6)
myWeight = weight(3)
mySecondWeight = weight(4)

mySupport = support((myPoint, mySecondPoint), (myWeight, mySecondWeight))

mySecondSupport = support((point(10,11), point(12,13)), (weight(5), weight(6)))

myMixedSubdivision = mixed_subdivision((mySupport, mySecondSupport))

export myMixedSubdivision

end
