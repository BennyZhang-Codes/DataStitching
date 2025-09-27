
function brain_hophantom3D(
    objbrain  ::BrainPhantom                 ; # PhantomType
    axis      ::String         = "axial"     , # orientation
    ss        ::Int64          = 4           , # undersample
    location  ::Vector{Int64}  = [160, 200]  , # [slice] for one slice, [start, end] for multiple slices
    b0_type   ::Symbol         = :fat        , # load B0 map
    csm_type  ::Symbol         = :real_32cha , # coil type
    csm_nCoil ::Int64          = 32          ,  
    verbose   ::Bool           = false       ,
) :: HO_Phantom

    class, size_phantom = load_phantom_mat(objbrain; axis=axis, ss=ss, location=location)

    if length(location) == 1 # 2D phantom
        nX = nY = nZ = 1
        if axis == "axial"
            nX, nY  = size(class)
        elseif axis == "coronal"
            nX, nZ  = size(class)
        elseif axis == "sagittal"
            nY, nZ  = size(class)
        end
    elseif length(location) == 2 # 3D phantom
        nX, nY, nZ  = size(class)
    end
    matrix_size = [nX, nY, nZ]
    spacing     = [objbrain.x*1e-3*ss, objbrain.y*1e-3*ss, objbrain.z*1e-3*ss];
    fov         = matrix_size .* spacing;
    Δx, Δy, Δz  = spacing

    # Generate 3D grid
    xyz = generate_grid(fov, matrix_size)
    x = xyz[:, :, :, 1]; 
    y = xyz[:, :, :, 2]; 
    z = xyz[:, :, :, 3]; 

    # Define spin property vectors
    T1, T2, T2s, ρ = SpinProperty_1p5T(class)

    # Define B0map vectors
    if b0_type == :fat
        Δw_fat = -220*2π
        Δw = (class.==93 )*Δw_fat .+ #FAT1
            (class.==209)*Δw_fat    #FAT2
    else    
        fieldmap = load_b0map(b0_type, size_phantom...; axis=axis, ss=ss, location=location, verbose=verbose)
        Δw = fieldmap*2π
    end

    # Define Coil-Sensitivity Map (CSM) vectors
    csm  = load_csm(csm_type, size_phantom..., csm_nCoil; axis=axis, ss=ss, location=location, verbose=verbose)
    nCha = size(csm)[end]

    # Define and return the Phantom struct
    x   = dropdims(  x; dims=Tuple(findall(size(  x) .== 1)));
    y   = dropdims(  y; dims=Tuple(findall(size(  y) .== 1)));
    z   = dropdims(  z; dims=Tuple(findall(size(  z) .== 1)));
    Δw  = dropdims( Δw; dims=Tuple(findall(size( Δw) .== 1)));
    csm = dropdims(csm; dims=Tuple(findall(size(csm) .== 1)));

    obj = HO_Phantom{Float64}(
        name = "brain3D_ss$(ss)_$(nX)x$(nY)x$(nZ)_Cha$(nCha)",
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
  