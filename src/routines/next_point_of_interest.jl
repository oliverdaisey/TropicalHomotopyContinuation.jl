include("../structs/mixed_cell_cone.jl")
include("ray_intersects_cone.jl")


"""
    next_point_of_interest(T::MixedCellTracker)

Given a mixed cell tracker, returns the next point of interest (either the next node, or the breaking point if it exists), along with the supports that change.

"""
function next_point_of_interest(T::MixedCellTracker)

    # make this code work with the old implementation (fix this after deadline)
    h = T.mixed_path
    pointer_index = 1
    fraction = QQFieldElem(0)

    # deal with the silly case that the mixed path only has one node left
    if pointer_index == length(h.pointers)
        return nothing, [] # this means that this mixed cell tracker is done
    end


    dual_path_pointers = h.pointers[pointer_index]
    n = length(dual_path_pointers)

    # get the lift at this time
    lift = vcat([lift_from_node_and_fraction(h.dualPaths[i], dual_path_pointers[i], fraction) for i in 1:n]...)

    println("lift = $(lift)")

    # get direction path is travelling in dual space
    next_dual_path_pointers = h.pointers[pointer_index+1]
    direction = vcat([h.dualPaths[i].nodes[next_dual_path_pointers[i]] for i in 1:n]...) ./ vcat([h.dualPaths[i].nodes[dual_path_pointers[i]] for i in 1:n]...)

    C_s = mixed_cell_cone(s)

    # starting at `lift` and moving in the direction `direction`, when do we hit facet of C_s?
    # we want to find the smallest t such that lift + t*direction is on the facet of C_s
    direction_data = [x.data for x in direction]
    t = ray_intersects_cone(mixed_cell_cone_to_polyhedron(C_s), lift, direction_data)


    if t > 1 || isnothing(t)
        # this means that the breaking point does not exist or takes us away from the path, so we should return the next node (no supports change)
        return h[2], []
    end

    # otherwise we are in the case that the breaking point exists, work out which supports change
    facet_point = lift + t*direction_data

    # get supports that change
    support_indices = []
    for i in 1:length(C_s.coefficients)
        if iszero(facet_point .* C_s.coefficients[i])
            # intersects this facet, supports that change are indices of nonzero entries of C_s.coefficients[i]
            push!(support_indices, findall(!iszero, C_s.coefficients[i]))
        end
    end

    # make sure support_indices is unique
    support_indices = unique(vcat(support_indices...))

    return TT.(facet_point), support_indices

end