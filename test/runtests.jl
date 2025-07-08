using Oscar
using TropicalHomotopyContinuation
using Test

@testset "TropicalHomotopyContinuation.jl" begin
    AbstractAlgebra.set_assertion_level(:TropicalHomotopyContinuation, 1)
    AbstractAlgebra.set_assertion_level(:TropicalHomotopyContinuationStart, 1)
    AbstractAlgebra.set_assertion_level(:TropicalHomotopyContinuationMove, 1)
    AbstractAlgebra.set_assertion_level(:TropicalHomotopyContinuationBergman, 1)
    AbstractAlgebra.set_assertion_level(:TropicalHomotopyContinuationJensen, 1)
    AbstractAlgebra.set_assertion_level(:TropicalHomotopyContinuationPerturb, 1)

    @testset "rigidity theory example" begin
        linearMatrix = matrix(QQ, [ 1  1  0  0  0  0  0  0  0  0;
                                    -1  0 -1 -1  0  0  0  0  0  0;
                                    0 -1  1  0 -1  0  0  0  0  0;
                                    0  0  0  1  1  0  0  0  0  0;
                                    0  0  0  0  0  1  1  0  0  0;
                                    0  0  0  0  0 -1  0 -1 -1  0;
                                    0  0  0  0  0  0 -1  1  0 -1;
                                    0  0  0  0  0  0  0  0  1  1])
        linearMatrix = echelon_form(linearMatrix, trim=true)
        M = matroid(linearMatrix)

        _,x,y = polynomial_ring(tropical_semiring(), "x"=>1:5, "y"=>1:5)
        F = vcat([xi*yi+0 for (xi, yi) in zip(x, y)],[x[1]+0])
        targetSupport = mixed_support(F)

        Oscar.randseed!(143)
        startingSupport, startingCell = starting_data(targetSupport, M)

        T = tracker(startingSupport, targetSupport, [startingCell], path=:straight_line)
        track!(T)
        Sigma = endgame(T)
        @test sum(Sigma[2]) == 4
    end

end
