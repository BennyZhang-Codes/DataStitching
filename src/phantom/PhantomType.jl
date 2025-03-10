abstract type PhantomType end



"""
    objb = BrainPhantom()
    objb = BrainPhantom(prefix::String; x::Float64=0.2, y::Float64=0.2, z::Float64=0.2)

The BrainPhantom. It contains parameters for selecting the Brain3D phantoms.

default: `brain3D_0.2_0.2_0.2.mat`.

# Arguments
- `prefix::String`: prefix of the Brain3D phantom.
- `x::Float64`:     x-voxelsize of the phantom.
- `y::Float64`:     y-voxelsize of the phantom.
- `z::Float64`:     z-voxelsize of the phantom.
- `name::String`:   name of the phantom, constructed from prefix, x, y, and z.
- `file::String`:   file name of the phantom, constructed from name.
- `matpath::String`: path of the phantom file `("*.mat")`.
# output
- `objb`: (`::BrainPhantom`) BrainPhantom struct.
"""
Base.@kwdef mutable struct BrainPhantom <: PhantomType 
    prefix::String = "brain3D"
    x::Float64 = 0.2
    y::Float64 = 0.2
    z::Float64 = 0.2
    name::String = "$(prefix)_$(z)_$(y)_$(x)"
    file::String = "$(name).mat"
    matpath::String = "$(@__DIR__)/mat/$(file)"
	BrainPhantom(prefix, x, y, z, name, file, matpath) = begin
		@assert isfile(matpath)  "The file $(file) does not exist."
		new(prefix, x, y, z, name, file, matpath)
	end
end

function BrainPhantom(prefix::String; x=0.2, y=0.2, z=0.2)
    return BrainPhantom(prefix=prefix, x=Float64(x), y=Float64(y), z=Float64(z))
end

Base.show(io::IO, b::BrainPhantom) = begin
	print(io, "BrainPhantom [$(b.prefix) | z=$(b.z) | y=$(b.y) | x=$(b.x)]")
end