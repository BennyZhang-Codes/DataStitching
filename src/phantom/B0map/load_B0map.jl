
"""
    B0map = load_B0map(B0_file::Symbol; axis="axial", ss=1, location=0.5)

Loads the B0map from the BrainPhantom package.

# Arguments
- `B0_file::Symbol`: The name of the B0map file to load. Must be one of `:B0` or `:B0_medianfiltered_r4`.
- `axis::String`: The axis along which to extract the B0map. Must be one of `"axial"`, `"coronal"`, or `"sagittal"`.
- `ss::Int`: The subsampling factor.
- `location::Float64`: The location along the specified axis at which to extract the B0map. Must be between 0 and 1.

# Returns
- `B0map::Array{Float64, 2}`: The B0map.
"""
function load_B0map(B0_file::Symbol; axis="axial", ss=1, location=0.5)
    @assert B0_file in [:B0, :B0_medianfiltered_r4] "B0_file must be one of the following: :B0, :B0_medianfiltered_r4"
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    if B0_file == :B0_medianfiltered_r4
        obj = BrainPhantom("brain3D_B0_medianfiltered_r4"; x=1, y=1, z=1)
    elseif B0_file == :B0
        obj = BrainPhantom("brain3D_B0"; x=1, y=1, z=1)
    end

    B0data = MAT.matread(obj.matpath)["data"];

    M, N, Z = size(B0data)
    if axis == "axial"
        loc   = Int32(ceil(Z*location))
        B0map = B0data[1:ss:end,1:ss:end, loc]
    elseif axis == "coronal"
        loc   = Int32(ceil(M*location))
        B0map = B0data[loc, 1:ss:end,1:ss:end]
    elseif axis == "sagittal"
        loc   = Int32(ceil(N*location))
        B0map = B0data[1:ss:end, loc,1:ss:end]
    end
	return B0map
end