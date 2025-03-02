module TropicalHomotopies

include("structs/point.jl")
include("structs/weight.jl")
include("structs/support.jl")
include("structs/mixed_support.jl")
include("structs/tracker.jl")
include("structs/cayley_embedding.jl")
include("structs/mixed_cell_cone.jl")

# Write your package code here.
p1 = point(0,0)
p2 = point(1,0)
p3 = point(0,1)
p4 = point(1,1)

p5=point(0,0)
p6=point(2,0)
p7=point(2,1)

mySupport = support([p1,p2,p3,p4], [weight(0), weight(0), weight(0), weight(0)])
mySecondSupport = support([p5,p6,p7], [weight(2), weight(1), weight(0)])

mixedSupport = mixed_support((mySupport, mySecondSupport))

cayley = cayley_embedding(mixedSupport)

candidate = mixed_support((support([p1,p2], [weight(0), weight(0)]), support([p6,p7], [weight(1), weight(0)])))

C = mixed_cell_cone(candidate, mixedSupport)

for facet in facets(C)
    println(facet)
end

println(candidate in mixed_cell_cone(candidate, mixedSupport))

polymakePolyhedron = convert(Polyhedron, mixed_cell_cone(candidate, mixedSupport))

display(Oscar.dim(polymakePolyhedron))

targetSupport = mixed_support((support([p1,p2,p3,p4], [weight(1), weight(0), weight(0), weight(0)]), support([p5,p6,p7], [weight(2), weight(1), weight(0)])))

T = tracker(mixedSupport, candidate, targetSupport)

pt, drift = tropical_intersection_point_and_drift(T)

println("pt = ", pt)

println("tropical drift = ", drift)

sop = mixed_subdivision(mixedSupport)

end