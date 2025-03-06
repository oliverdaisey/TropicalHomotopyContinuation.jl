module TropicalHomotopies

using Oscar

include("structs/point.jl")
include("structs/height.jl")
include("structs/support.jl")
include("structs/realisable_matroid.jl")
include("structs/chain_of_flats.jl")
include("structs/mixed_support.jl")
include("structs/mixed_cell.jl")
include("structs/tracker.jl")
include("structs/cayley_embedding.jl")
include("structs/mixed_cell_cone.jl")
include("jensen_move.jl")
include("bergman_move.jl")

p1 = point(0,0,0,2)
p2 = point(0,0,0,0)

p3=point(-1,-1,-1,0)
p4=point(0,0,0,0)

f3 = support([p1,p2],[0,0])
f4 = support([p3,p4],[0,2])

f3Target = support([p1,p2],[0,0])
f4Target = support([p3,p4],[0,-1])

mixedSupport = mixed_support((f3,f4))
targetSupport = mixed_support((f3Target,f4Target))

M = matroid(Oscar.matrix(QQ, [-1 -2 -3 -4; 11 13 15 19]))

chainOfFlats = chain_of_flats(M, [[4]])
println("loopless face = ", loopless_face(chainOfFlats))

candidate = mixed_cell(mixed_support((support([p1,p2],[0,0]),support([p3,p4],[0,3]))), chainOfFlats)

T = tracker(mixedSupport, [candidate], [targetSupport])

println("Computing mixed cell cone")
C = mixed_cell_cone(candidate, mixedSupport)
println("Mixed cell cone = ", C)

# show the heights of T before the bergman move
println(" ")
println("heights before Bergman move")
for p in points(ambient_support(T))
    println("p = ", p, " in support ", findfirst(S -> p in S, supports(ambient_support(T))), " height = ", ambient_support(T)[p])
end
println(" ")

println("Mixed cells before Bergman move")
println(" ")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println("loopless face = ", loopless_face(chain_of_flats(σ)))
    println(" ")
end

println("* bergman move *")
bergman_move!(T)

# show the heights of T after the bergman move
println(" ")
println("heights after Bergman move")
for p in points(ambient_support(T))
    println("p = ", p, " in support ", findfirst(S -> p in S, supports(ambient_support(T))), " height = ", ambient_support(T)[p])
end
println(" ")

println("Mixed cells after Bergman move")
println(" ")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println("loopless face = ", loopless_face(chain_of_flats(σ)))
    println(" ")
end


end