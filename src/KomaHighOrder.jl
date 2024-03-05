module KomaHighOrder


import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty, Base.copy
# IMPORT PACKAGES
using CUDA, Interpolations, PlotlyJS
using ThreadsX
# Printing
using ProgressMeter
using Reexport
@reexport using KomaMRI

import KomaMRI.KomaMRICore: SimulationMethod, SpinStateRepresentation
@reexport import MAT
import Functors: @functor

include("datatypes/Sequence.jl")
include("mrd/signal_to_raw_data.jl")

include("simulation/DiscreteSequence.jl") # include HO_Sequence
include("simulation/Simulate.jl")
include("simulation/SimulatorCore.jl")
include("simulation/Bloch/BlochSimulatitonMethod.jl")
include("simulation/Bloch/HO1BlochSimulatitonMethod.jl")
include("plot/plot_seq.jl")
include("plot/plot_seqd.jl")

export HO_signal_to_raw_data
export HO_DiscreteSequence, HO_discretize
export HO_Sequence
export HO_simulate
export HO_plot_seq
export HO_plot_seqd, HO_plot_hoseqd

#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaHighOrder
