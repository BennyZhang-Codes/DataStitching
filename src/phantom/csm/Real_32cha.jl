"""
    csm = csm_Real_32cha(nX, nY; verbose=false)

# Description
    This function loads the Real Coil-Sensitivity Map (32 channels) and resizes it to the desired size (nX, nY, 32).

# Arguments
- `nX::Int64`: number of voxels in the x-direction
- `nY::Int64`: number of voxels in the y-direction
- `verbose::Bool=false`: whether to print progress messages

# Returns
- `smap::Array{ComplexF64, 3}`: the coil-sensitivity map

```julia
>>> smap = csm_Real_32cha(200, 200)
```
"""
function csm_Real_32cha(nX::Int64, nY::Int64; verbose::Bool=false)
    if verbose
        @info "loading Real Coil-Sensitivity Map (32 channels)..."
    end
    sensitivity = MAT.matread("$(@__DIR__)/coilsensmap_32cha.mat")["coilsensmap_32cha"]
    out = imresize(sensitivity, (nX, nY, 32))
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
    out = out./ norm
    out[isnan.(out)] .= 0 + 0im
    return out
end

"""
    csm = csm_Real_32cha(nX, nY, nZ; verbose=false)

# Description
    This function loads the Real Coil-Sensitivity Map (32 channels) and resizes it to the desired size (nX, nY, nZ, 32).

# Arguments
- `nX::Int64`: number of voxels in the x-direction
- `nY::Int64`: number of voxels in the y-direction
- `nZ::Int64`: number of voxels in the z-direction
- `verbose::Bool=false`: whether to print progress messages

# Returns
- `smap::Array{ComplexF64, 4}`: the coil-sensitivity map

```julia
>>> smap = csm_Real_32cha(200, 200, 70)
```
"""
function csm_Real_32cha(nX::Int64, nY::Int64, nZ::Int64; verbose::Bool=false)
    if verbose
        @info "loading Real Coil-Sensitivity Map (32 channels)..."
    end
    sensitivity = MAT.matread("$(@__DIR__)/coilsensmap_32cha_3d.mat")["csm"]
    out = imresize(sensitivity, (nX, nY, nZ, 32))
    norm = sqrt.(sum(abs.(out) .^ 2, dims=4))
    out = out./ norm
    out[isnan.(out)] .= 0 + 0im
    return out
end