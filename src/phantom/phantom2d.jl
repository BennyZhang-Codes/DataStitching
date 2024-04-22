import KomaMRI.KomaMRIBase: brain_phantom2D
# import Base.show

Base.@kwdef mutable struct brain2D <: PhantomType 
    name::String = "0.2 mm isotropic brain phantom"
    file::String = "brain3D_0.2.mat"
end
export brain2D

"""
    obj = brain_phantom2D(p::PhantomType; axis="axial", ss=4, location=0.5, loadB0map=true, mask_idx=1, Nparts=1)
Creates a two-dimensional brain Phantom struct.

# References
- B. Aubert-Broche, D.L. Collins, A.C. Evans: "A new improved version of the realistic
    digital brain phantom" NeuroImage, in review - 2006
- B. Aubert-Broche, M. Griffin, G.B. Pike, A.C. Evans and D.L. Collins: "20 new digital
    brain phantoms for creation of validation image data bases" IEEE TMI, in review - 2006
- https://brainweb.bic.mni.mcgill.ca/brainweb

# Keywords
- `axis`: (`::String`, `="axial"`, opts=[`"axial"`, `"coronal"`, `"sagittal"`]) orientation of the phantom
- `ss`: (`::Integer`, `=4`) subsampling parameter in all axis
- `location`: (`::Float64`, `=0.5`) location of the phantom in the Z-axis
- `loadB0map`: (`::Bool`, `=true`) load B0 map
- `mask_idx`: (`::Integer`, `=1`) mask index
- `Nparts`: (`::Integer`, `=1`) number of partitions (split in fan shape)
```
"""
function brain_phantom2D(
    p::PhantomType;         # phantom type
    axis::String="axial",   # orientation
    ss::Int64=4,            # undersample
    location::Float64=0.5,  # relative location in the Z-axis
    loadB0map::Bool=true,   # load B0 map
    mask_idx::Int64=1,      # mask index
    Nparts::Int64=1,        # number of partitions (split in fan shape)
    ) :: Phantom
    path = (@__DIR__) * phantom_dict[:path]
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert isfile(path*"/$(p.file)") "File $(p.file) does not exist in $(path)"
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    @assert 1 <= mask_idx <= Nparts "mask_idx must be between 1 and Nparts"

    data = MAT.matread(path*"/$(p.file)")["data"]
    
    M, N, Z = size(data)
    if axis == "axial"
        loc = Int32(ceil(Z*location))
        class = data[1:ss:end,1:ss:end, loc]
    elseif axis == "coronal"
        loc = Int32(ceil(M*location))
        class = data[loc, 1:ss:end,1:ss:end]   
    elseif axis == "sagittal"
        loc = Int32(ceil(N*location))
        class = data[1:ss:end, loc,1:ss:end]
    end

    # Define spin position vectors
    Δx = .2e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

    # Define spin property vectors
    T2 = (class.==23)*329 .+ #CSF
        (class.==46)*83 .+ #GM
        (class.==70)*70 .+ #WM
        (class.==93)*70 .+ #FAT1
        (class.==116)*47 .+ #MUSCLE
        (class.==139)*329 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*70 .+ #FAT2
        (class.==232)*329 .+ #DURA
        (class.==255)*70 #MARROW
    T2s = (class.==23)*58 .+ #CSF
        (class.==46)*69 .+ #GM
        (class.==70)*61 .+ #WM
        (class.==93)*58 .+ #FAT1
        (class.==116)*30 .+ #MUSCLE
        (class.==139)*58 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*61 .+ #FAT2
        (class.==232)*58 .+ #DURA
        (class.==255)*61 .+#MARROW
        (class.==255)*70 #MARROW
    T1 = (class.==23)*2569 .+ #CSF
        (class.==46)*833 .+ #GM
        (class.==70)*500 .+ #WM
        (class.==93)*350 .+ #FAT1
        (class.==116)*900 .+ #MUSCLE
        (class.==139)*569 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*500 .+ #FAT2
        (class.==232)*2569 .+ #DURA
        (class.==255)*500 #MARROW
    ρ = (class.==23)*1 .+ #CSF
        (class.==46)*.86 .+ #GM
        (class.==70)*.77 .+ #WM
        (class.==93)*1 .+ #FAT1
        (class.==116)*1 .+ #MUSCLE
        (class.==139)*.7 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*.77 .+ #FAT2
        (class.==232)*1 .+ #DURA
        (class.==255)*.77 #MARROW
    if loadB0map    
        B0map = brain_phantom2D_B0map(; axis=axis, ss=1, location=location)
        B0map = imresize(B0map, size(class))
        Δw = B0map*2π
    else
        Δw_fat = -220*2π
        Δw = (class.==93)*Δw_fat .+ #FAT1
            (class.==209)*Δw_fat    #FAT2
    end

	T1 = T1*1e-3
	T2 = T2*1e-3
	T2s = T2s*1e-3

    mask = get_fan_mask(M, N, Nparts)

    # Define and return the Phantom struct
    obj = Phantom{Float64}(
        name = "brain2D_$(axis)_ss$(ss)_location$(location)-$(loc)_$(mask_idx)/$(Nparts)",
		x = y[ρ.!=0 .&& mask.==mask_idx],
		y = x[ρ.!=0 .&& mask.==mask_idx],
		z = 0*x[ρ.!=0 .&& mask.==mask_idx],
		ρ = ρ[ρ.!=0 .&& mask.==mask_idx],
		T1 = T1[ρ.!=0 .&& mask.==mask_idx],
		T2 = T2[ρ.!=0 .&& mask.==mask_idx],
		T2s = T2s[ρ.!=0 .&& mask.==mask_idx],
		Δw = Δw[ρ.!=0 .&& mask.==mask_idx],
    )
	return obj
end
  
function brain_phantom2D_B0map(; axis="axial", ss=1, location=0.5)
    path = (@__DIR__) * phantom_dict[:path]
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert isfile(path*"/$(phantom_dict[:brain2d_B0])") "file is not found in $(path)"
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    B0data = MAT.matread(path*"/$(phantom_dict[:brain2d_B0])")["data"];

    M, N, Z = size(B0data)
    if axis == "axial"
        loc = Int32(ceil(Z*location))
        B0map = B0data[1:ss:end,1:ss:end, loc]
    elseif axis == "coronal"
        loc = Int32(ceil(M*location))
        B0map = B0data[loc, 1:ss:end,1:ss:end]
    elseif axis == "sagittal"
        loc = Int32(ceil(N*location))
        B0map = B0data[1:ss:end, loc,1:ss:end]
    end
	return B0map
end
