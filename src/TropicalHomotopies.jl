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
include("homotopy.jl")
include("starting_system.jl")


# M = matroid(Oscar.matrix(QQ, [1 1 0 0 0 1 1 0 0 0; -1 0 -1 -1 0 -1 0 -1 -1 0; 0 -1 1 0 -1 0 -1 1 0 -1; 0 0 0 1 1 0 0 0 1 1]))
# targetSupports = Support[]
# # M = matroid(Oscar.matrix(QQ, [1 1 0 0 0 0 0 0 0 ; -1 0 -1 -1 0 0 0 0 0 ; 0 -1 1 0 -1 0 0 0 0 ; 0 0 0 1 1 0 0 0 0 ; 0 0 0 0 0 1 1 0 0 ; 0 0 0 0 0 -1 0 -1 -1 ; 0 0 0 0 0 0 -1 1 0 ; 0 0 0 0 0 0 0 0 1 1]))
# for i in [1,2,3,4,5]
#         pi = point([n in [i,i+5] ? 1 : 0 for n in 1:10])
#         p0 = point([0 for n in 1:10])
#         push!(targetSupports, support([pi,p0],[0,0]))
# end

# targetSupport = mixed_support(targetSupports)

# for S in supports(targetSupport)
#     println(S)
# end
# # # NEW EXAMPLE

# Δ, σ = starting_data(targetSupport, M)

# T = tracker(Δ, targetSupport, [σ])

# display(stable_intersection(T))

end 