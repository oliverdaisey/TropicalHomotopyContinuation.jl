module TropicalHomotopyContinuation

using Oscar

# the following functions in Oscar will be extended:
import Oscar:
    dot,
    flats,
    matrix,
    stable_intersection

include("structs/point.jl")
include("structs/height.jl")
include("structs/support.jl")
include("structs/realisable_matroid.jl")
include("structs/flat.jl")
include("structs/chain_of_flats.jl")
include("structs/chain_of_flats_cone.jl")
include("structs/mixed_support.jl")
include("structs/mixed_cell.jl")
include("structs/logger.jl")
include("structs/tracker.jl")
include("structs/cayley_embedding.jl")
include("structs/mixed_cell_cone.jl")
include("jensen_flip.jl")
include("bergman_flip.jl")
include("homotopy.jl")
include("starting_system.jl")
include("exports.jl")
include("interface.jl")

function __init__()
    # scope for overall homotopies
    add_verbosity_scope(:TropicalHomotopyContinuation)
    add_assertion_scope(:TropicalHomotopyContinuation)

    # scope for starting data
    add_verbosity_scope(:TropicalHomotopyContinuationStart)
    add_assertion_scope(:TropicalHomotopyContinuationStart)

    # scope for the general move
    add_verbosity_scope(:TropicalHomotopyContinuationMove)
    add_assertion_scope(:TropicalHomotopyContinuationMove)

    # scope for Bergman related functions (time and move)
    add_verbosity_scope(:TropicalHomotopyContinuationBergman)
    add_assertion_scope(:TropicalHomotopyContinuationBergman)

    # scope for Jensen related functions (time and move)
    add_verbosity_scope(:TropicalHomotopyContinuationJensen)
    add_assertion_scope(:TropicalHomotopyContinuationJensen)

    # scope for perturbation
    add_verbosity_scope(:TropicalHomotopyContinuationPerturb)
    add_assertion_scope(:TropicalHomotopyContinuationPerturb)
end

end
