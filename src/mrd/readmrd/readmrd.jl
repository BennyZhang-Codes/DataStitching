#=
In order to handle the rawdata in ISMRMRD format, we write some functions to extract the necessary information from the header and data.
The following functions/datas are included:
    - `mrddims`: a list of the dimensions of the MRD data
    - `get_kinfo`: get the fov, matrix size and k-space shape
    - `get_ksize`: get the shape of the k-space data
    - `get_kdata`: get the k-space data
=#

"""
a list of the dimensions of the MRD data:
- `nCha`: channel, Number of channels
- `nZ`: kspace_encode_step_2, Partition encoding
- `nY`: kspace_encode_step_1, Phase encoding line
- `nX`: readout, number of readout points 
- `nAvg`: average, Signal average
- `nSli`: slice, Slice number (multi-slice 2D)
- `nCon`: contrast, Echo number in multi-echo
- `nPha`: phase, Cardiac phase
- `nRep`: repetition, Counter in repeated/dynamic acquisitions
- `nSet`: set, Sets of different preparation, e.g. flow encoding, diffusion weighting
- `nSeg`: segment, Counter for segmented acquisitions
"""
mrddims = ["nCha", "nZ", "nY", "nX", "nAvg", "nSli", "nCon", "nPha", "nRep", "nSet", "nSeg"]
export mrddims

include("get_size.jl")
export get_isize, get_ksize

include("get_data.jl")
export get_kdata