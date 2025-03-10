include("SignalOp.jl")

using LinearOperators
abstract type HOOp{T} <: AbstractLinearOperator{T} end

include("HighOrderOp.jl")
include("HighOrderOp_v1.jl")
include("HighOrderOp_v1p1.jl")
include("HighOrderOp_v2.jl")
include("HighOrderOp_v3.jl")