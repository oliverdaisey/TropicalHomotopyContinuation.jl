using Revise
using Oscar
using TropicalHomotopyContinuation
include("examples/rigidity_theory_example.jl")


M = matrix(QQ,[1 0 1; 0 1 1])
Mmat = matroid(RealisableMatroid(M,2))
F = TropicalHomotopyContinuation.flat(Mmat, Set{Int}([1]))
G = TropicalHomotopyContinuation.flat(matroid(RealisableMatroid(M,2)), Set{Int}([1]))
G in [F]



# unrefined reduced flat: 6,8,9,10
# refined reduced flats: 6,8 and 9,10
solution = QQFieldElem[-344, -635, -635, -1145//2, -635, -604, -897//2, -604, -1209//2, -1209//2]




w = QQFieldElem[-344, -635, -635, -572, -635, -604, -434, -604, -604, -590]
u = QQFieldElem[0, 0, 0, 0, 0, 0, -1//2, 0, 0, -1//2]
tBergman = 28

jensenStarLinearEquationMatrix = matrix(QQ, [2 0 0 0 0 -2 0 0 0 0; 0 0 0 0 0 0 2 0 0 0; 0 0 -2 0 0 0 0 2 0 0; 0 0 0 -2 0 0 0 0 2 0; 0 0 0 0 0 0 0 0 0 -2; 1 0 0 0 0 0 0 0 0 0])
jensenTrail = convex_hull([w + tBergman * u], [u], transpose(Oscar.kernel(jensenStarLinearEquationMatrix, side=:right)))

bergmanConeRayColumnMatrix = matrix(QQ,[1 1 1 1 1 1; 0 0 0 0 0 1; 0 0 0 0 0 1; 0 0 1 1 1 1; 0 0 0 0 0 1; 0 0 0 0 1 1; 0 1 1 1 1 1; 0 0 0 0 1 1; 0 0 0 0 1 1; 0 0 0 1 1 1])
chainOfFlatsCone = positive_hull(transpose(bergmanConeRayColumnMatrix), ones_matrix(QQ, 1,nrows(bergmanConeRayColumnMatrix)))
Oscar.dim(intersect(jensenTrail, polyhedron(chainOfFlatsCone)))==0 # true, so far so good



affineLinearEquationsLHS = matrix(QQ,[0 -1 1 0 0 0 0 0 0 0; 0 -1 0 0 1 0 0 0 0 0; 0 0 0 0 0 -1 0 1 0 0; 0 0 0 0 0 -1 0 0 1 0; 2 0 0 0 0 -2 0 0 0 0; 0 0 0 0 0 0 2 0 0 0; 0 0 -2 0 0 0 0 2 0 0; 0 0 0 -2 0 0 0 0 2 0; 0 0 0 0 0 0 0 0 0 -2; 1 0 0 0 0 0 0 0 0 0])
affineLinearEquationsRHS = QQFieldElem[0, 0, 0, 0, 520, -897, 62, -64, 1209, -344]
canSolve, solution, kernelGenerators = Oscar.can_solve_with_solution_and_kernel(affineLinearEquationsLHS, affineLinearEquationsRHS; side=:right)
canSolve == true # true
solution == QQFieldElem[-344, -635, -635, -572, -635, -604, -897//2, -604, -604, -1209//2]
solution in jensenTrail # true, so far so good


A = affineLinearEquationsLHS[1:4,:]
b = affineLinearEquationsRHS[1:4]

sigmaBerg = polyhedron((zero_matrix(QQ,0,10),QQFieldElem[]),affine_hull(polyhedron(chainOfFlatsCone)))
deltaBerg = polyhedron((A,b),(A,b))
sigmaBerg == deltaBerg # true
solution in sigmaBerg # true, so far so good


tropicalPoint = (10000000*(w+tBergman*u)+solution)/10000001
sortperm(tropicalPoint,rev=true)

# Unrefined chain:
# ∅ ⊊ {1} ⊊ {1, 7} ⊊ {1, 4, 7} ⊊ {1, 4, 6, 7, 8, 9, 10} ⊊ {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
# Old flipped chains:
# ∅ ⊊ {1} ⊊ {1, 7} ⊊ {1, 4, 7} ⊊ {1, 4, 6, 7, 8} ⊊ {1, 4, 6, 7, 8, 9, 10} ⊊ {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
# New flipped chains:
# ∅ ⊊ {1} ⊊ {1, 7} ⊊ {1, 4, 7} ⊊ {1, 4, 6, 7, 8} ⊊ {1, 4, 6, 7, 8, 9, 10} ⊊ {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
# ∅ ⊊ {1} ⊊ {1, 7} ⊊ {1, 4, 7} ⊊ {1, 4, 7, 10} ⊊ {1, 4, 6, 7, 8, 9, 10} ⊊ {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} # redundant, but coincides with matrix used to construct chainOfFlatsCone


(1000000*(w+tBergman*u)+solution)/1000001 in chainOfFlatsCone # false, so w+tBergman*u+epsilon*solution not in chainOfFlatsCone

raysUnrefinedChain = transpose(bergmanConeRayColumnMatrix[:,[1,2,3,5,6]])
extraRay = bergmanConeRayColumnMatrix[:,4]
bergmanConeStar = positive_hull([extraRay], raysUnrefinedChain)
(1000000*(w+tBergman*u)+solution)/1000001 in bergmanConeStar
!(solution in bergmanConeStar)
w+tBergman*u in bergmanConeStar
