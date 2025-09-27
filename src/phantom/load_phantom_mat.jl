function load_phantom_mat(
    objbrain  :: BrainPhantom                ;  # PhantomType
    axis      :: String         = "axial"    ,  # orientation
    ss        :: Int64          = 4          ,  # undersample
    location  :: Vector{Int64}  = [160, 200] ,  # [slice] for one slice, [start, end] for multiple slices
)
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert length(location) in (1,2) "location must have length 1 (single slice) or 2 (slice range)"

    data = MAT.matread(objbrain.matpath)["data"]

    nY, nX, nZ = size(data)
    size_phantom = (nX, nY, nZ)

    if length(location) == 1 # 2D phantom
        idx = location[1]
        if axis == "axial"
            @assert 0 < idx <= nZ "slice index out of range"
            class = data[1:ss:end, 1:ss:end, idx]
        elseif axis == "coronal"
            @assert 0 < idx <= nY "slice index out of range"
            class = data[idx, 1:ss:end, 1:ss:end]   
        elseif axis == "sagittal"
            @assert 0 < idx <= nX "slice index out of range"
            class = data[1:ss:end, idx, 1:ss:end]
        end
     elseif length(location) == 2 # 3D phantom
        start_idx, end_idx = location
        @assert start_idx <= end_idx "start index must be <= end index"
        if axis == "axial"
            @assert 0 < start_idx <= nZ && 0 < end_idx <= nZ "slice range out of bounds"
            class = data[1:ss:end, 1:ss:end, start_idx:ss:end_idx]
        elseif axis == "coronal"
            @assert 0 < start_idx <= nY && 0 < end_idx <= nY "slice range out of bounds"
            class = data[start_idx:ss:end_idx, 1:ss:end, 1:ss:end]   
        elseif axis == "sagittal"
            @assert 0 < start_idx <= nX && 0 < end_idx <= nX "slice range out of bounds"
            class = data[1:ss:end, start_idx:ss:end_idx, 1:ss:end]
        end
    else
        throw(ArgumentError("location must have length 1 (single slice) or 2 (slice range)"))
    end
    return class, size_phantom
end