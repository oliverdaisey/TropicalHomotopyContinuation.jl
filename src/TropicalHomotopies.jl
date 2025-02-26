module TropicalHomotopies

using Oscar
include("structs/point.jl")
include("structs/weight.jl")
include("structs/support.jl")

# Write your package code here.
myPoint = point(3,4)
mySecondPoint = point(5,6)
myWeight = weight(3)
mySecondWeight = weight(4)

mySupport = support((myPoint, mySecondPoint), (myWeight, mySecondWeight))
display(mySupport)
end
