"""
    shape = get_ksize(raw::RawAcquisitionData)

# Description
    get the k-space size of the raw acquisition data.

# Arguments
- `raw::RawAcquisitionData`: the raw acquisition data.

# Returns
- `shape::Tuple`: [nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg]

# Example
```julia-repl
julia> raw = RawAcquisitionData(ISMRMRDFile("path/to/file.mrd"))
julia> shape = get_ksize(raw)
julia> nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape
```
"""
function get_ksize(raw::RawAcquisitionData)
    profiles = raw.profiles
    nCha = nZ = nY = nX = nAvg = nSli = nCon = nPha = nRep = nSet = nSeg = 1 
    for idx in eachindex(profiles)
        profile = profiles[idx]
        header = profile.head
        # Stuff into the buffer
        # kspace_encode_step_1, kspace_encode_step_2, average, slice, contrast, phase, repetition, set, segment
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
        # data[:, Par, Lin, :, Avg, Sli, Con, Pha, Rep, Set, Seg] = transpose(profile.data)
        nCha = max(nCha, Cha)
        nZ   = max(nZ  ,   Z)
        nY   = max(nY  ,   Y)
        nX   = max(nX  ,   X)
        nAvg = max(nAvg, Avg)
        nSli = max(nSli, Sli)
        nCon = max(nCon, Con)
        nPha = max(nPha, Pha)
        nRep = max(nRep, Rep)
        nSet = max(nSet, Set)
        nSeg = max(nSeg, Seg)
        # println((Cha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg))
    end
    shape = nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg
    return shape
end

"""
    fov, matrix, shape = get_kinfo(raw::RawAcquisitionData)

# Description
    get the fov, matrix size, k-space shape of the raw acquisition data.

# Arguments
- `raw::RawAcquisitionData`: the raw acquisition data.

# Returns
- `fov::Tuple`: (eFOVx, eFOVy, eFOVz, rFOVx, rFOVy, rFOVz)
- `matrix::Tuple`: (eNx, eNy, eNz, rNx, rNy, rNz)
- `shape::Tuple`: (nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg)
"""
function get_kinfo(raw::RawAcquisitionData)
    Flag_IgnoreSeg = true
    header = raw.params

    # Matrix size
    eNx, eNy, eNz = raw.params["encodedSize"]
    rNx, rNy, rNz = raw.params["reconSize"]

    # Field of View
    eFOVx, eFOVy, eFOVz = raw.params["encodedFOV"]
    rFOVx, rFOVy, rFOVz = raw.params["reconFOV"]

    # Number of Slices, Reps, Contrasts, etc.
    # kspace_encode_step_1, kspace_encode_step_2, average, slice, contrast, phase, repetition, set, segment, user
    nCha = raw.params["receiverChannels"]
    nAvg = 1 + header["enc_lim_average"].maximum
    nSli = 1 + header["enc_lim_slice"].maximum
    nCon = 1 + header["enc_lim_contrast"].maximum
    nPha = 1 + header["enc_lim_phase"].maximum
    nRep = 1 + header["enc_lim_repetition"].maximum
    nSet = 1 + header["enc_lim_set"].maximum
    nSeg = Flag_IgnoreSeg ?  1 : header["enc_lim_segment"].maximum + 1

    fov    = (eFOVx, eFOVy, eFOVz, rFOVx, rFOVy, rFOVz)
    matrix = (eNx, eNy, eNz, rNx, rNy, rNz)
    shape  = (nCha, rNz, max(eNy, rNy), rNx, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg)
    return fov, matrix, shape
end

