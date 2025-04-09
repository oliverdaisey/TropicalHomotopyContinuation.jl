mutable struct Logger

    maxMixedCells::Int
    numberOfMoves::Int
    numberOfJensenMoves::Int
    numberOfBergmanMoves::Int
    numberOfSimultaneousBergmanAndJensenMoves::Int
    numberOfDivergedMixedCells::Int
    numberOfRebases::Int

end

function logger()
    return Logger(0, 0, 0, 0, 0, 0, 0)
end

function Base.show(io::IO, L::Logger)

    print(io, "Max mixed cells $(L.maxMixedCells), number of moves $(L.numberOfMoves), number of jensen moves $(L.numberOfJensenMoves), number of bergman moves $(L.numberOfBergmanMoves), number of simultaneous Bergman and Jensen moves: $(L.numberOfSimultaneousBergmanAndJensenMoves), number of rebases: $(L.numberOfRebases), number of diverged mixed cells $(L.numberOfDivergedMixedCells)")
end