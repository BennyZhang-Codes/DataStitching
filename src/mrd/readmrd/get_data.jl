"""
    kdata = get_kdata(raw::RawAcquisitionData)

# Description
    get the k-space data from the raw acquisition data.

# Arguments
- `raw::RawAcquisitionData`: the raw acquisition data.

# Returns
- `kdata::Array{Complex{Float32}, 11}`: k-space data with dimensions [nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg]

# Example
```julia-repl
julia> raw = RawAcquisitionData(ISMRMRDFile("path/to/file.mrd"))
julia> kdata = get_kdata(raw)
```
"""
function get_kdata(raw::RawAcquisitionData, shape::Tuple)
    profiles = raw.profiles
    kdata = Array{Complex{Float32},length(shape)}(undef,shape);
    for idx in eachindex(profiles)
        profile = profiles[idx]
        header = profile.head
        # Stuff into the buffer
        Cha = Int(profile.head.active_channels)
        Z   = 1 + header.idx.kspace_encode_step_2
        Y   = 1 + header.idx.kspace_encode_step_1
        X   = Int(profile.head.number_of_samples)
        Avg = 1 + header.idx.average
        Sli = 1 + header.idx.slice
        Con = 1 + header.idx.contrast
        Pha = 1 + header.idx.phase
        Rep = 1 + header.idx.repetition
        Set = 1 + header.idx.set
        Seg = 1 + header.idx.segment
        kdata[:, Z, Y, :, Avg, Sli, Con, Pha, Rep, Set, Seg] = transpose(profile.data);
    end
    return kdata
end

