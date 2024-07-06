import KomaMRI.KomaMRIBase: brain_phantom2D
# import Base.show

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
    objbrain::BrainPhantom;        # PhantomType
    axis::String="axial",          # orientation
    ss::Int64=4,                   # undersample
    location::Float64=0.5,         # relative location in the Z-axis
    B0_type::Symbol=:fat,          # load B0 map
    B0_file::Symbol=:B0,
    maxOffresonance::Float64=125.,
    coil_type::Symbol=:fan,        # coil type
    coil_idx::Int64 =1,            # coil index
    Nparts::Int64   =1,            # number of partitions (split in fan shape)
    Npartsx::Int64  =1,            # number of partitions in the X-axis
    Npartsy::Int64  =1,            # number of partitions in the Y-axis
    overlap::Real   =1,            # overlap between fan coils
    relative_radius::Real=1.5,     # relative radius of the coil
    ) :: Phantom
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert B0_type in [:real, :fat, :quadratic] "B0_type must be one of the following: :real, :fat, :quadratic"
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    @assert coil_type in [:fan, :rect, :birdcage, :real] "coil_type must be one of the following: :fan, :rect, :birdcage, :real"
    if coil_type == :rect
        Nparts = Npartsx * Npartsy
    elseif coil_type == :real
        Nparts = 32
    end
    @assert 1 <= coil_idx <= Nparts "coil_idx must be between 1 and Nparts"
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

    # Define spin position vectors
    Δx = .2e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

    # Define spin property vectors
    T2 = (class.==23 )*329 .+ #CSF
         (class.==46 )*83  .+ #GM
         (class.==70 )*70  .+ #WM
         (class.==93 )*70  .+ #FAT1
         (class.==116)*47  .+ #MUSCLE
         (class.==139)*329 .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*70  .+ #FAT2
         (class.==232)*329 .+ #DURA
         (class.==255)*70     #MARROW
    T2s =(class.==23 )*58  .+ #CSF
         (class.==46 )*69  .+ #GM
         (class.==70 )*61  .+ #WM
         (class.==93 )*58  .+ #FAT1
         (class.==116)*30  .+ #MUSCLE
         (class.==139)*58  .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*61  .+ #FAT2
         (class.==232)*58  .+ #DURA
         (class.==255)*61  .+ #MARROW
         (class.==255)*70     #MARROW
    T1 = (class.==23 )*2569 .+ #CSF
         (class.==46 )*833 .+ #GM
         (class.==70 )*500 .+ #WM
         (class.==93 )*350 .+ #FAT1
         (class.==116)*900 .+ #MUSCLE
         (class.==139)*569 .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*500 .+ #FAT2
         (class.==232)*2569 .+ #DURA
         (class.==255)*500    #MARROW
    ρ =  (class.==23 )*1   .+ #CSF
         (class.==46 )*.86 .+ #GM
         (class.==70 )*.77 .+ #WM
         (class.==93 )*1   .+ #FAT1
         (class.==116)*1   .+ #MUSCLE
         (class.==139)*.7  .+ #SKIN/MUSCLE
         (class.==162)*0   .+ #SKULL
         (class.==185)*0   .+ #VESSELS
         (class.==209)*.77 .+ #FAT2
         (class.==232)*1   .+ #DURA
         (class.==255)*.77    #MARROW

    T1  = T1  * 1e-3
    T2  = T2  * 1e-3
    T2s = T2s * 1e-3

    if B0_type == :real    
        fieldmap = brain_phantom2D_B0map(B0_file; axis=axis, ss=1, location=location)
        fieldmap = imresize(fieldmap, size(class))
        Δw = fieldmap*2π
    elseif B0_type == :fat
        Δw_fat = -220*2π
        Δw = (class.==93 )*Δw_fat .+ #FAT1
             (class.==209)*Δw_fat    #FAT2
    elseif B0_type == :quadratic
        fieldmap = quadraticFieldmap(size(class)...,maxOffresonance)[:,:,1]
        Δw = fieldmap*2π
    end


    if coil_type == :fan
        smap = get_fan_mask(M, N, Nparts; overlap=overlap)[:,:,coil_idx]
    elseif coil_type == :rect
        smap = get_rect_mask(M, N, Npartsx, Npartsy)[:,:,coil_idx]
    elseif coil_type == :birdcage
        smap = BirdcageSensitivity(M, N, Nparts; relative_radius=Float64(relative_radius))[:,:,coil_idx]
    elseif coil_type == :real
        smap = RealCoilSensitivity_32cha(M, N)[:,:,coil_idx]
    end

    # Define and return the Phantom struct
    obj = Phantom{Float64}(
        name = "brain2D_$(axis)_ss$(ss)_$(M)x$(N)_location$(location)-$(loc)_$(coil_idx)/$(Nparts)",
		x   =    y[ρ.!=0],
		y   =    x[ρ.!=0],
		z   =  0*x[ρ.!=0],
		ρ   =    ρ[ρ.!=0],
		T1  =   T1[ρ.!=0],
		T2  =   T2[ρ.!=0],
		T2s =  T2s[ρ.!=0],
		Δw  =   Δw[ρ.!=0],
        Dλ1 = real(smap)[ρ.!=0],
        Dλ2 = imag(smap)[ρ.!=0],
        Dθ  = abs.(smap)[ρ.!=0],
    )
	return obj
end
  
function brain_phantom2D_B0map(B0_file::Symbol; axis="axial", ss=1, location=0.5)
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
