module KomaHighOrder


import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty, Base.copy
# IMPORT PACKAGES
using CUDA, Interpolations, PlotlyJS
using Reexport
@reexport using KomaMRI
@reexport import MAT

export greet
greet() = print("Hello World!")

include("mrd/HO_signal_to_raw_data.jl")
include("simulation/HO_simulate.jl")
include("datatypes/HO_Sequence.jl")

include("simulation/HO_DiscreteSequence.jl") # include HO_Sequence
include("ui/DisplayFunctions.jl")

export HO_signal_to_raw_data
export HO_DiscreteSequence, HO_discretize
export HO_Sequence
export HO_simulate
export HO_plot_seq, HO_plot_seqd

#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaHighOrder
