module KomaHighOrder


import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty, Base.copy, Base.show
# IMPORT PACKAGES
using CUDA, Interpolations, PlotlyJS
using ThreadsX
# Printing
using ProgressMeter
using Reexport
using MRIBase
@reexport using KomaMRI

import KomaMRI.KomaMRICore: SimulationMethod, SpinStateRepresentation

@reexport import MAT
import Functors: @functor

include("datatypes/Sequence.jl")
include("phantom/phantom2d.jl")
include("phantom/phantom3d.jl")

include("simulation/DiscreteSequence.jl") # include HO_Sequence
include("simulation/Simulate.jl")
include("simulation/SimulatorCore.jl")
include("simulation/Bloch/BlochSimulatitonMethod.jl")
include("simulation/Bloch/BlochHighOrderSimulatitonMethod.jl")
include("plot/plot_seq.jl")
include("plot/plot_seqd.jl")
include("plot/plot_kspace.jl")
include("reconstruction/recon_2d.jl")

include("mrd/signal_to_raw_data.jl")
include("mrd/HO_signal_to_raw_data.jl")
# include("mrd/AcquisitionData.jl")

export HO_DiscreteSequence
export HO_Sequence
export simulate
export plot_seq
export plot_seqd, plot_hoseqd
# export signal_to_raw_data

export reconstruct_2d_image

export info
#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module KomaHighOrder
