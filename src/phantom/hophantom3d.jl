
function brain_hophantom3D(
    objbrain::BrainPhantom         ;  # PhantomType
    ss::Int64           = 4        ,  # undersample
    start_end           = [160, 200],

    db0_type::Symbol    = :fat     ,  # load B0 map
    db0_file::Symbol    = :B0      ,  # determines the *.mat file of the B0 map

    csm_type::Symbol    = :real_32cha ,  # coil type
    csm_nCoil::Int64    = 32        ,  
    verbose::Bool       = false     ,
) :: HO_Phantom
    @assert db0_type in [:real, :fat, :quadratic] "db0_type must be one of the following: :real, :fat, :quadratic"

    class, _ = load_phantom_mat(objbrain; ss=ss, start_end=start_end)

    M, N, Z = size(class)
    matrix_size = [M, N, Z]
    spacing     = [objbrain.x*1e-3*ss, objbrain.y*1e-3*ss, objbrain.z*1e-3*ss];
    fov         = matrix_size .* spacing;
    Δx, Δy, Δz = spacing

    # Generate 3D grid
    xyz = generate_grid(fov, matrix_size)
    x = xyz[:, :, :, 1];
    y = xyz[:, :, :, 2];
    z = xyz[:, :, :, 3];

    # Define spin property vectors
    T1, T2, T2s, ρ = SpinProperty_1p5T(class)

    # Define B0map vectors
    if db0_type == :real    
        fieldmap = load_B0map(db0_file; ss=1, start_end=start_end)
        fieldmap = imresize(fieldmap, size(class))
        Δw = fieldmap*2π
    elseif db0_type == :fat
        Δw_fat = -220*2π
        Δw = (class.==93 )*Δw_fat .+ #FAT1
            (class.==209)*Δw_fat    #FAT2
    end

    # Define Coil-Sensitivity Map (CSM) vectors
    csm = load_csm(csm_type, M, N, Z, csm_nCoil; verbose=verbose)
    csm_nCoil = size(csm)[end]

    # Define and return the Phantom struct
    obj = HO_Phantom{Float64}(
        name = "brain3D_ss$(ss)_$(M)x$(N)x$(Z)_Cha$(csm_nCoil)",
        x   =    y[ρ.!=0],
        y   =    x[ρ.!=0],
        z   =    z[ρ.!=0],
        ρ   =    ρ[ρ.!=0],
        T1  =   T1[ρ.!=0],
        T2  =   T2[ρ.!=0],
        T2s =  T2s[ρ.!=0],
        Δw  =   Δw[ρ.!=0],
        csm =  csm[ρ.!=0,:],
    )
	return obj
end
  