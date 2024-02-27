module PolyPlanning

using LinearAlgebra
using OSQP
using SparseArrays
using Infiltrator
using GLMakie

include("solvers.jl")
include("poly_functions.jl")

export ConvexPolygon2D

end # module PolyPlanning
