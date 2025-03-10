


convert_ifft(x;dims=[1,2]) = fftshift(ifft(ifftshift(x,dims),dims),dims)*prod(size(x)[dims])
convert_fft(x;dims=[1,2])  = fftshift( fft(ifftshift(x,dims),dims),dims)/prod(size(x)[dims])


"""
    get_kdata(raw::RawAcquisitionData)
This function takes a `RawAcquisitionData` object and returns the k-space data as a 3D array.
# Arguments
`raw` :: `RawAcquisitionData`
# Returns
`k` :: `Array{ComplexF64, 3}`, (Nx, Ny, Ncoils)
# Example
```julia-repl
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


function recon_2d_ifft(raw::RawAcquisitionData)
    k = get_kdata(raw)
    return convert_ifft(k)
end