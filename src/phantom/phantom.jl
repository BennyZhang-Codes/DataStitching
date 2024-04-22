abstract type PhantomType end

function info(s::Phantom)
	print("Phantom[name = $(s.name) | spins = $(length(s.x)) | x = $(minimum(s.x)*1e2):$(maximum(s.x)*1e2) cm | y = $(minimum(s.y)*1e2):$(maximum(s.y)*1e2) cm | z = $(minimum(s.z)*1e2):$(maximum(s.z)*1e2) cm ]")
    print("\n")
end

Base.show(io::IO, s::PhantomType) = begin
    print(io, "name = $(s.name)\n")
    print(io, "file = $(s.file)\n")
end

phantom_dict = Dict{Symbol, String}(
    :path        => "/mat",
    :brain2d     => "brain3D_0.2.mat",
    :brain2d_B0  => "brain3D_B0map_1.0.mat",
    :brain3d_171 => "(171, 191)_brain3D_(400, 362, 434)_(0.025, 0.5, 0.5).mat",
    :brain3d_285 => "(285, 305)_brain3D_(400, 362, 434)_(0.025, 0.5, 0.5).mat",
)


export info
export phantom_dict

include("phantom2d.jl")
include("phantom3d.jl")
include("PhantomReference.jl")

export PhantomReference, brain_phantom2D_B0map