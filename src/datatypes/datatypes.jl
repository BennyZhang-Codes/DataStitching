#=
This file contains the definitions of the data types used in the package.
In HighOrderMRI.jl, we extended the three data types based on KomaMRI.jl:
    1. HO_Sequence: sequence with field dynamics.
    2. HO_DiscreteSequence: discretized HO_Sequence.
    3. HO_Phantom: Phantom with Coil-Sensitivity maps.
=#

include("HO_Sequence.jl")         
include("HO_DiscreteSequence.jl") 
include("HO_Phantom.jl")


export HO_DiscreteSequence
export HO_Sequence
export HO_Phantom
