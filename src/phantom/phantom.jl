# definition of the abstract PhantomType and BrainPhantom struct
include("PhantomType.jl")
export BrainPhantom

include("SpinProperty.jl")
export SpinProperty_1p5T

export load_phantom_mat
function load_phantom_mat(
    objbrain::BrainPhantom;        # PhantomType
    axis::String="axial",          # orientation
    ss::Int64=4,                   # undersample
    location::Float64=0.5,         # relative location in the slice direction
)
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    data = MAT.matread(objbrain.matpath)["data"]
    M, N, Z = size(data)
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
    return class, loc
end


# HO_Phantom type, support coil-sensitivity map (CSM)
include("hohantom2d.jl")
export brain_hophantom2D

# include the Coil-Sensitivity Map (CSM) module
include("csm/csm.jl")

include("B0map/B0map.jl")



# functions to generate different types of phantoms
include("phantom2d.jl")
include("phantom2d_reference.jl")
include("phantom3d.jl")

export brain_phantom2D_reference

# function to print the information of a Phantom object
function info(s::Phantom)
	print("Phantom[$(s.name) | nSpin=$(length(s.x)) | x=$(minimum(s.x)*1e2):$(maximum(s.x)*1e2) cm | y=$(minimum(s.y)*1e2):$(maximum(s.y)*1e2) cm | z=$(minimum(s.z)*1e2):$(maximum(s.z)*1e2) cm ]")
    print("\n")
end

export info

