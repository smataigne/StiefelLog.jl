module Manifold

using LinearAlgebra, SkewLinearAlgebra, MatrixEquations
import LinearAlgebra as LA
export 
    Orthonormal,
    StiefelVector,
    TangentVector,
    randorthonormal,
    orthogonalize!,
    orthogonalize,
    isorthonormal,
    gramschmidt,
    gramschmidt!,
    buildcomplement,
    addcomplement,
    dims,
    projection,
    projection!,
    gettangent


include("orthonormal.jl")
include("stiefel.jl")
include("exp.jl")
include("log.jl")
include("pathstraightening.jl")

end

