import KomaMRI.KomaMRIBase: brain_phantom2D
# import Base.show

"""
    obj = brain_phantom2D(p::BrainPhantom; axis="axial", ss=5, location=0.5, B0_type=:fat, B0_file=:B0, maxOffresonance=125., csmtype=:fan, coil_idx=1, nCoil=1, overlap=1, relative_radius=1.5)
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
    ss::Int64=5,                   # undersample
    location::Float64=0.5,         # relative location in the Z-axis

    B0_type::Symbol=:fat,          # load B0 map
    B0_file::Symbol=:B0,
    maxOffresonance::Float64=125.,

    csmtype::Symbol=:fan,          # coil type
    coil_idx::Int64 =1,            # coil index
    nCoil::Int64   =1,             # number of Coils
    overlap::Real   =1,            # overlap between fan coils, for csm_Fan_binary
    relative_radius::Real=1.5,     # relative radius of the coil, for csm_Birdcage
    ) :: Phantom
    @assert B0_type in [:real, :fat, :quadratic] "B0_type must be one of the following: :real, :fat, :quadratic"
    @assert 1 <= coil_idx <= nCoil "coil_idx must be between 1 and $(nCoil)"

    class = load_phantom_mat(objbrain; axis=axis, ss=ss, location=location)
    # Define spin position vectors
    Δx = .2e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

    # Define spin property vectors
    T1, T2, T2s, ρ = SpinProperty_1p5T(class)

    # Define B0map vectors
    if B0_type == :real    
        fieldmap = load_B0map(B0_file; axis=axis, ss=1, location=location)
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
    
    # Define Coil-Sensitivity Map (CSM) vectors
    csm = load_csm(csmtype, M, N, nCoil; overlap=overlap, relative_radius=relative_radius)
    _,_,nCoil = size(csm)
    csm = csm[:,:,coil_idx]

    # Define and return the Phantom struct
    obj = Phantom{Float64}(
        name = "brain2D_$(axis)_ss$(ss)_$(M)x$(N)_location$(location)_$(coil_idx)/$(nCoil)",
		x   =    y[ρ.!=0],
		y   =    x[ρ.!=0],
		z   =  0*x[ρ.!=0],
		ρ   =    ρ[ρ.!=0],
		T1  =   T1[ρ.!=0],
		T2  =   T2[ρ.!=0],
		T2s =  T2s[ρ.!=0],
		Δw  =   Δw[ρ.!=0],
        Dλ1 = real(csm)[ρ.!=0],
        Dλ2 = imag(csm)[ρ.!=0],
        Dθ  = abs.(csm)[ρ.!=0],
    )
	return obj
end
  

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
    return class
end
