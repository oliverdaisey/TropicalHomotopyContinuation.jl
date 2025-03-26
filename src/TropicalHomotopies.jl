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

# M = matroid(Oscar.matrix(QQ, [-1 -2 -3 -4; 11 13 15 19]))

# p1 = point(0,0,0,1)
# p2 = point(0,0,0,0)

# p3 = point(0,0,0,0)
# p4 = point(0,1,1,0)
# p5 = point(0,2,0,0)
# p6 = point(0,0,2,0)
# p7 = point(0,2,2,0)

# t = -3
# hts = [0, -3 + t, 0 + 2*t, -4, -4 + 2*t]
# println("hts = ", hts)
# f1 = support([p1,p2],[0,0])
# f2 = support([p3,p4,p5,p6,p7],hts)
# t = 3
# hts = [0, -3 + t, 0 + 2*t, -4, -4 + 2*t]
# f2Target = support([p3,p4,p5,p6,p7],hts)
# targetSupport = mixed_support((f1,f2Target))

# w = starting_data(targetSupport, M)

# println("w = ", w)

include("../examples/bergman_move_example.jl")

end 