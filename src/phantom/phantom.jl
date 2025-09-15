# definition of the abstract PhantomType and BrainPhantom struct
include("PhantomType.jl")
export BrainPhantom

include("SpinProperty.jl")
export SpinProperty_1p5T

include("load_phantom_mat.jl")
export load_phantom_mat

# HO_Phantom type, support coil-sensitivity map (CSM)
include("hophantom2d.jl")
export brain_hophantom2D

# include the BrainPhantom module
include("hophantom3d.jl")
export brain_hophantom3D

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

