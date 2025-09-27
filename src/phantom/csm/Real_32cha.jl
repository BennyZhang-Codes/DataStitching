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
    csm = csm_Real_32cha(nX, nY, nZ; axis="axial", ss=4, location=[160,200], verbose=false)

# Description
    This function loads the Real Coil-Sensitivity Map (32 channels) and resizes it to the desired size (nX, nY, nZ, 32).

# Arguments
- `nX::Int64`: number of voxels of the phantom in the x-direction  
- `nY::Int64`: number of voxels of the phantom in the y-direction  
- `nZ::Int64`: number of voxels of the phantom in the z-direction 
- `axis::String="axial"`: orientation of slices, must be one of `"axial"`, `"coronal"`, `"sagittal"`  
- `ss::Int64=4`: undersampling factor for spatial dimensions  
- `location::Vector{Int64}=[160,200]`:  
    - if length is 1 → index of a single slice  
    - if length is 2 → `[start, end]` slice range  
- `verbose::Bool=false`: whether to print progress messages  

# Returns
- `smap::Array{ComplexF64, 4}`: the coil-sensitivity map

```julia
>>> smap = csm_Real_32cha(200, 200, 70; axis="axial", ss=4, location=[180])
```
"""
function csm_Real_32cha(
    nX       :: Int64                       , 
    nY       :: Int64                       , 
    nZ       :: Int64                       ; 
    axis     :: String         = "axial"    ,  # orientation
    ss       :: Int64          = 4          ,  # undersample
    location :: Vector{Int64}  = [160, 200] ,  # [slice] for one slice, [start, end] for multiple slices
    verbose  :: Bool           = false      ,
)
    @assert length(location) in (1,2) "location must have length 1 (single slice) or 2 (slice range)"

    if verbose
        @info "loading Real Coil-Sensitivity Map (32 channels)..."
    end
    csm = MAT.matread("$(@__DIR__)/csm_4p0.mat")["csm"]
    y, x, z, nCha = size(csm)

    res = 4e-3; # [m], resolution of the original csm
    res_phantom = 0.5e-3; # [m], resolution of the phantom

    range_y = range(0.5, y - 0.5, length=54) * res
    range_x = range(0.5, x - 0.5, length=45) * res
    range_z = range(0.5, z - 0.5, length=45) * res

    if length(location) == 1 # 2D phantom
        idx = location[1]
        if axis == "axial"
            @assert 0 < idx <= nZ "slice index out of range"
            x_new = res_phantom * (collect(1:ss:nX) .- 0.5)
            y_new = res_phantom * (collect(1:ss:nY) .- 0.5)
            z_new = res_phantom * [idx]
        elseif axis == "coronal"
            @assert 0 < idx <= nY "slice index out of range"
            x_new = res_phantom * (collect(1:ss:nX) .- 0.5)
            y_new = res_phantom * [idx]
            z_new = res_phantom * (collect(1:ss:nZ) .- 0.5)
        elseif axis == "sagittal"
            @assert 0 < idx <= nX "slice index out of range"
            x_new = res_phantom * [idx]
            y_new = res_phantom * (collect(1:ss:nY) .- 0.5)
            z_new = res_phantom * (collect(1:ss:nZ) .- 0.5)
        end
    elseif length(location) == 2 # 3D phantom
        start_idx, end_idx = location
        @assert start_idx <= end_idx "start index must be <= end index"
        if axis == "axial"
            @assert 0 < start_idx <= nZ && 0 < end_idx <= nZ "slice range out of bounds"
            x_new = res_phantom * (collect(1:ss:nX) .- 0.5)
            y_new = res_phantom * (collect(1:ss:nY) .- 0.5)
            z_new = res_phantom * (collect(start_idx:ss:end_idx) .- 0.5)
        elseif axis == "coronal"
            @assert 0 < start_idx <= nY && 0 < end_idx <= nY "slice range out of bounds"
            x_new = res_phantom * (collect(1:ss:nX) .- 0.5)
            y_new = res_phantom * (collect(start_idx:ss:end_idx) .- 0.5)
            z_new = res_phantom * (collect(1:ss:nZ) .- 0.5)
        elseif axis == "sagittal"
            @assert 0 < start_idx <= nX && 0 < end_idx <= nX "slice range out of bounds"
            x_new = res_phantom * (collect(start_idx:ss:end_idx) .- 0.5)
            y_new = res_phantom * (collect(1:ss:nY) .- 0.5)
            z_new = res_phantom * (collect(1:ss:nZ) .- 0.5)
        end
    else
        throw(ArgumentError("location must have length 1 (single slice) or 2 (slice range)"))
    end

    itp = interpolate(csm, BSpline(Linear()))
    itp = Interpolations.scale(itp, range_y, range_x, range_z, 1:nCha)
    eitp = extrapolate(itp, Flat())
    csm_out = [eitp(y,x,z,cha) for y in y_new, x in x_new, z in z_new, cha in 1:nCha]

    norm = sqrt.(sum(abs.(csm_out) .^ 2, dims=4))
    csm_out = csm_out./ norm
    csm_out[isnan.(csm_out)] .= 0 + 0im
    return csm_out
end