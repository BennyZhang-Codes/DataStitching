"""
    HighOrderMRI.jl

A Julia package for high-order MRI simulation and reconstruction.

1. For simulation: we extend the KomaMRI.jl to realize high-order and parallel MRI simulation.
2. For reconstruction: we extend the MRIReco.jl to realize high-order MRI reconstruction with off-resonance reconstruction.
3. Synchronization: we also implemented a model-based delay estimation method for synchronization between the MRI data and the trajectory data.

"""
module HighOrderMRI

import Base.*, Base.+, Base.-, Base./, Base.vcat, Base.size, Base.abs, Base.getproperty, Base.copy, Base.show
# IMPORT PACKAGES


using PyPlot#, PlotlyJS  

using Parameters
using ThreadsX

using ProgressMeter
using Reexport

using CUDA, Interpolations, Statistics

using AbstractNFFTs, NFFTTools
using FFTW: fftshift, ifftshift, fft, ifft
import ImageTransformations: imresize

@reexport import MAT

using MRIReco, RegularizedLeastSquares
import MRISimulation: quadraticFieldmap   

@reexport using KomaMRI
import KomaMRI.KomaMRICore: SimulationMethod, SpinStateRepresentation
import Functors: @functor


include("utils/utils.jl")
include("datatypes/datatypes.jl")
include("phantom/phantom.jl")
include("simulation/simulation.jl")
include("mrd/mrd.jl")
include("reconstruction/recon.jl")
include("load/load.jl")
include("synchronization/synchronization.jl")

include("plot/plot.jl")

#Package version
using Pkg
__VERSION__ = VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module HighOrderMRI
