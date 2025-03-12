using HighOrderMRI
using Test
using MRISimulation

@info "testing... InterpTraj.jl"
include("test_InterpTraj.jl")

@info "testing... csm.jl"
include("test_csm.jl")

@info "testing... HighOrderOp.jl"
include("test_HighOrderOp.jl")

@info "testing... HOPhantom.jl"
include("test_HOPhantom.jl")