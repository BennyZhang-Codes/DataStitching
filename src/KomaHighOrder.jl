module KomaHighOrder


import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty, Base.copy, Base.show
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
include("phantom/phantom.jl")

include("simulation/DiscreteSequence.jl") # include HO_Sequence
include("simulation/Simulate.jl")
include("simulation/SimulatorCore.jl")
include("simulation/Bloch/BlochSimulatitonMethod.jl")
include("simulation/Bloch/BlochHighOrderSimulatitonMethod.jl")
include("plot/plot_seq.jl")
include("plot/plot_seqd.jl")


export HO_DiscreteSequence
export HO_Sequence
export simulate
export plot_seq
export plot_seqd, plot_hoseqd
export HO_signal_to_raw_data

#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaHighOrder
