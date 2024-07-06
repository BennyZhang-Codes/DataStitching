
function brain_hophantom2D(
    objbrain::BrainPhantom;        # PhantomType
    axis::String="axial",          # orientation
    ss::Int64=5,                   # undersample
    location::Float64=0.5,         # relative location in the Z-axis

    B0_type::Symbol=:fat,          # load B0 map
    B0_file::Symbol=:B0,
    maxOffresonance::Float64=125.,

    csmtype::Symbol=:fan,          # coil type
    nCoil::Int64   =1,            # number of partitions (split in fan shape)
    # overlap::Real   =1,            # overlap between fan coils
    # relative_radius::Real=1.5,     # relative radius of the coil
    ) :: HO_Phantom
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert B0_type in [:real, :fat, :quadratic] "B0_type must be one of the following: :real, :fat, :quadratic"
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    
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

    csm = load_csm(csmtype, M, N, nCoil)
    _,_,nCoil = size(csm)

    # Define and return the Phantom struct
    obj = HO_Phantom{Float64}(
        name = "brain2D_$(axis)_ss$(ss)_$(M)x$(N)_location$(location)-$(loc)_Cha$(nCoil)",
		x   =    y[ρ.!=0],
		y   =    x[ρ.!=0],
		z   =  0*x[ρ.!=0],
		ρ   =    ρ[ρ.!=0],
		T1  =   T1[ρ.!=0],
		T2  =   T2[ρ.!=0],
		T2s =  T2s[ρ.!=0],
		Δw  =   Δw[ρ.!=0],
        csm =  csm[ρ.!=0,:],
    )
	return obj
end
  