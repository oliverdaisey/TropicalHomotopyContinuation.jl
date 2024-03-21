# This file includes all the files in the src directory in order of dependency.

include("structs/dual_type.jl")
include("structs/support.jl")
include("structs/tropical_pluecker_vector.jl")
include("structs/dual_cell.jl")
include("structs/dual_path.jl")
include("structs/mixed_cell.jl")
include("structs/mixed_path.jl")
include("structs/mixed_cell_tracker.jl")
include("routines/cayley_embedding.jl")
include("structs/mixed_cell_cone.jl")
include("routines/ray_intersects_cone.jl")
include("routines/subdivision_of_points.jl")
include("routines/tropical_drift.jl")
include("routines/stable_intersection_point.jl")
include("routines/perturb_dual_vector.jl")
include("routines/next_point_of_interest.jl")
include("routines/is_dual_cell_candidate.jl")
include("routines/inflation_deflation.jl")
include("routines/dual_subdivisions.jl")