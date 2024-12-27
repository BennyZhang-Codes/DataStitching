module KomaHighOrder


import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty, Base.copy, Base.show
# IMPORT PACKAGES
using CUDA, Interpolations, Statistics
import ImageTransformations: imresize
import MRISimulation: quadraticFieldmap   
using PyPlot, PlotlyJS  

using Parameters
using ThreadsX

using ProgressMeter
using Reexport

using MRIReco, RegularizedLeastSquares
@reexport using KomaMRI

using FFTW: fftshift, ifftshift, fft, ifft
import KomaMRI.KomaMRICore: SimulationMethod, SpinStateRepresentation

@reexport import MAT
import Functors: @functor

using AbstractNFFTs
using NFFTTools

include("utils/utils.jl")
include("datatypes/datatypes.jl")
include("phantom/phantom.jl")
include("simulation/simulation.jl")
include("plot/plot.jl")
include("mrd/mrd.jl")
include("reconstruction/recon.jl")
include("example/example.jl")
include("synchronization/synchronization.jl")

#Package version, KomaMRI.__VERSION__
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])


end # module KomaHighOrder
