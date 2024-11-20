
function brain_hophantom2D(
    objbrain::BrainPhantom         ;  # PhantomType
    axis::String        = "axial"  ,  # orientation
    ss::Int64           = 5        ,  # undersample
    location::Float64   = 0.5      ,  # relative location in the Z-axis

    db0_type::Symbol    = :fat     ,  # load B0 map
    db0_file::Symbol    = :B0      ,  # determines the *.mat file of the B0 map
    db0_max::Float64    = 125.     ,  # for quadraticFieldmap

    csm_type::Symbol    = :fan     ,  # coil type
    csm_nCoil::Int64    = 1        ,  # number of partitions (split in fan shape)
    csm_nRow            = nothing  ,
    csm_nCol            = nothing  ,
    csm_nBlock          = 3        ,
    csm_overlap::Real   = 0        ,  # overlap between fan coils, for csm_Fan_binary
    csm_radius::Real    = 1.5      ,  # relative radius of the coil, for csm_Birdcage

    verbose::Bool       = false    ,
) :: HO_Phantom
    @assert db0_type in [:real, :fat, :quadratic] "db0_type must be one of the following: :real, :fat, :quadratic"

    class, loc = load_phantom_mat(objbrain; axis=axis, ss=ss, location=location)
    # Define spin position vectors
    Δx = Δy = objbrain.x*1e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δy #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δy:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

    # Define spin property vectors
    T1, T2, T2s, ρ = SpinProperty_1p5T(class)

    # Define B0map vectors
    if db0_type == :real    
        fieldmap = load_B0map(db0_file; axis=axis, ss=1, location=location)
        fieldmap = imresize(fieldmap, size(class))
        Δw = fieldmap*2π
    elseif db0_type == :fat
        Δw_fat = -220*2π
        Δw = (class.==93 )*Δw_fat .+ #FAT1
             (class.==209)*Δw_fat    #FAT2
    elseif db0_type == :quadratic
        fieldmap = quadraticFieldmap(size(class)...,db0_max)[:,:,1]
        Δw = fieldmap*2π
    end

    # Define Coil-Sensitivity Map (CSM) vectors
    csm = load_csm(csm_type, M, N, csm_nCoil; overlap=csm_overlap, 
        relative_radius=csm_radius, nRow=csm_nRow, nCol=csm_nCol, nBlock=csm_nBlock, verbose=verbose)
    _,_,csm_nCoil = size(csm)

    # Define and return the Phantom struct
    obj = HO_Phantom{Float64}(
        name = "brain2D_$(axis)_ss$(ss)_$(M)x$(N)_loc$(location)-$(loc)_Cha$(csm_nCoil)",
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
  