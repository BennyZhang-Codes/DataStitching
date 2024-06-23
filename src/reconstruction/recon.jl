include("EncodingOperators/SignalOp.jl")
include("EncodingOperators/HighOrderOp.jl")
include("recon_2d.jl")
include("reconstrution.jl")
include("FFT.jl")

export convert_fft, convert_ifft
export reconstruct_2d_image

export get_kdata, reconstruct_2d_ifft


"""
    get_kdata(raw::RawAcquisitionData)
This function takes a `RawAcquisitionData` object and returns the k-space data as a 3D array.
# Arguments
`raw` :: `RawAcquisitionData`
# Returns
`k` :: `Array{ComplexF64, 3}`, (Nx, Ny, Ncoils)
# Example
```julia
>>> k = get_kdata(raw)
>>> plot_img(abs.(k[:,:,1]).^0.1)
>>> plot_img(abs.(reconstruct_ifft(k))[:,:,1])
````
"""
function get_kdata(raw::RawAcquisitionData)
    numProfiles = length(raw.profiles)
    numSamplesPerProfile, numChannels = size(raw.profiles[1].data)
    k = Array{ComplexF64, 3}(undef, numSamplesPerProfile, numProfiles, numChannels)

    for i in 1:numProfiles
        k[:, i, :] = raw.profiles[i].data
    end
    return k
end


function reconstruct_2d_ifft(raw::RawAcquisitionData)
    k = get_kdata(raw)
    return convert_ifft(k)
end