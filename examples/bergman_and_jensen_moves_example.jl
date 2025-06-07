using TropicalHomotopyContinuation
import Oscar.uniform_matroid # do not import all of Oscar to avoid name clashes

M = matroid(Oscar.matrix(QQ, [-1 -2 -3 -4; 11 13 15 19]))

p1 = TropicalHomotopyContinuation.point(0,0,0,1)
p2 = TropicalHomotopyContinuation.point(0,0,0,0)

p3 = TropicalHomotopyContinuation.point(0,0,0,0)
p4 = TropicalHomotopyContinuation.point(0,1,1,0)
p5 = TropicalHomotopyContinuation.point(0,2,0,0)
p6 = TropicalHomotopyContinuation.point(0,0,2,0)
p7 = TropicalHomotopyContinuation.point(0,2,2,0)

t = -3
hts = [0, -3 + t, 0 + 2*t, -4, -4 + 2*t]
println("hts = ", hts)
f1 = TropicalHomotopyContinuation.support([p1,p2],[0,0])
f2 = TropicalHomotopyContinuation.support([p3,p4,p5,p6,p7],hts)



t = 3
hts = [0, -3 + t, 0 + 2*t, -4, -4 + 2*t]
f2Target = TropicalHomotopyContinuation.support([p3,p4,p5,p6,p7],hts)

mixedSupport = mixed_support((f1,f2))
targetSupport = mixed_support((f1,f2Target))

# candidate one corresponds to the jensen moves
f2Active = TropicalHomotopyContinuation.support([p5, p7], [6, 2])
chainOfFlats = chain_of_flats(M, [[3]])
candidateOne = mixed_cell(mixed_support((f1, f2Active)), chainOfFlats)

# candidate two corresponds to the bergman moves
f2Active = TropicalHomotopyContinuation.support([p6, p7], [-4, -10])
chainOfFlats = chain_of_flats(M, [[2]])
candidateTwo = mixed_cell(mixed_support((f1, f2Active)), chainOfFlats)

# check that the intersection points are correct
T = tracker(mixedSupport, targetSupport, [candidateOne, candidateTwo], path=:straight_line)

# @time display(stable_intersection(T))

println("candidate one")
w, u = tropical_intersection_point_and_drift(T, candidateOne)

println("candidate two")
w, u = tropical_intersection_point_and_drift(T, candidateTwo)

println("-----------------------------------------------")
println("Printing mixed cells")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println(" ")
end

move!(T)


println("-----------------------------------------------")
println("Printing new mixed cells after move 1")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println(" ")
end

move!(T)

println("-----------------------------------------------")
println("Printing new mixed cells after move 2")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println(" ")
end

move!(T)

println("-----------------------------------------------")
println("Printing new mixed cells after move 3")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println(" ")
end

move!(T)

println("-----------------------------------------------")
println("Printing new mixed cells after move 4")
for σ in mixed_cells(T)
    println(σ)
    println("pt and drift = ", tropical_intersection_point_and_drift(T, σ))
    println(" ")
end