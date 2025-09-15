function load_phantom_mat(
    objbrain::BrainPhantom;        # PhantomType
    axis::String="axial",          # orientation
    ss::Int64=4,                   # undersample
    start_end = [160, 200],        # only used for 3D phantom
    location::Float64=0.5,         # relative location in the slice direction
)
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    data = MAT.matread(objbrain.matpath)["data"]
    M, N, Z = size(data)
    if Z == 1 # 2D phantom
        if axis == "axial"
            loc   = Int32(ceil(Z*location))
            class = data[1:ss:end,1:ss:end, loc]
        elseif axis == "coronal"
            loc   = Int32(ceil(M*location))
            class = data[loc, 1:ss:end,1:ss:end]   
        elseif axis == "sagittal"
            loc   = Int32(ceil(N*location))
            class = data[1:ss:end, loc,1:ss:end]
        end
    else # 3D phantom
        loc = Int32(0)
        class = data[1:ss:end, 1:ss:end, start_end[1]:ss:start_end[2]]
    end
    return class, loc
end