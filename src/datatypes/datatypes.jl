#=
This file contains the definitions of the data types used in the package.
In KomaHighOrder.jl, we built the three data types based on KomaMRI.jl:
    1. HO_Sequence: sequence with dynamic fields.
    2. HO_DiscreteSequence: discretized HO_Sequence.
    3. HO_Phantom: phantom with Coil-Sensitivity maps.
=#

include("HO_Sequence.jl")         
include("HO_DiscreteSequence.jl") 
include("HO_Phantom.jl")
include("SphericalHarmonics.jl")
# grid for reconstruction
include("grid/grid.jl")

export HO_DiscreteSequence
export HO_Sequence
export HO_Phantom
export SphericalHarmonics