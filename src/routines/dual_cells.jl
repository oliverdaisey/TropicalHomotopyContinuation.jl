using Oscar

"""
    dual_cells(S::Support{<:DualType}, c::Vector{Oscar.TropicalSemiringElem})

Create a vector of all the maximal dual cells from a support (of exponent vectors) and a lift.

TODO: Return *all* dual cells, not just the maximal ones.
"""
function dual_cells(S::Support{<:DualType}, c::Vector{Oscar.TropicalSemiringElem})


    dualSubdivision = subdivision_of_points_workaround(S.points, c)
    polyhedralComplex = polyhedral_complex(Oscar.pm_object(dualSubdivision).POLYHEDRAL_COMPLEX)

    dualCells = []
    for i in 1:dim(polyhedralComplex)
        for p in polyhedra_of_dim(polyhedralComplex, i)
            Vector{Int}.(collect(vertices(p)))
        end
    end

    

    return [DualCell(S, active_support, c) for active_support in maximal_cells(dual_subdivision)]
end