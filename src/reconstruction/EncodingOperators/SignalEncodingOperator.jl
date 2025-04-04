include("SignalOp.jl")

using LinearOperators
abstract type HOOp{T} <: AbstractLinearOperator{T} end

# the formal version HighOrderOp, the extended signal model with field dynamics and off-resonance
include("HighOrderOp.jl")

# developing versions
include("HighOrderOp_v1.jl")
include("HighOrderOp_v1p1.jl")
include("HighOrderOp_v2.jl")
include("HighOrderOp_v3.jl")