

function brain_phantom2D_reference(
    objbrain::BrainPhantom        ,
    key::Symbol                   ,
    target_fov::Tuple{T,T}        ,  # [mm] target FOV
    target_res::Tuple{T,T}        ;  # [mm] target resolution
    axis              = "axial"   , 
    ss::Int64         = 3         , 
    location::T       = 0.8       , 

    db0_type::Symbol  = :fat      ,  # load B0 map
    db0_file::Symbol  = :B0       ,  # determines the *.mat file of the B0 map
    db0_max::T        = 125.      ,  # for quadraticFieldmap

    csm_type::Symbol  = :fan      ,  # coil type
    csm_nCoil::Int64  = 1         ,  # number of partitions (split in fan shape)
    csm_nRow          = nothing   ,
    csm_nCol          = nothing   ,
    csm_nBlock        = 3         ,
    csm_overlap::Real = 0         ,  # overlap between fan coils, for csm_Fan_binary
    csm_radius::Real  = 1.5       ,  # relative radius of the coil, for csm_Birdcage

    verbose::Bool     = false     ,
) where {T<:Float64}



    @assert key in [:ρ, :T2, :T2s, :T1, :raw, :headmask, :brainmask, :Δw, :csm] "key must be ρ, T2, T2s, T1, Δw, raw, headmask or brainmask"
    @assert db0_type in [:real, :fat, :quadratic] "db0_type must be one of the following: :real, :fat, :quadratic"
    @assert 0 <= location <= 1 "location must be between 0 and 1"

    Δx, Δy = [objbrain.x, objbrain.y] .* ss;   # phantom resoultion, original resolution multiplied by undersampling factor ss.
    fov_x, fov_y = target_fov
    res_x, res_y = target_res
    center_range = (Int64(ceil(fov_x / (Δx))), Int64(ceil(fov_y / (Δy)))) 
    target_size = (Int64(ceil(fov_x / res_x)), Int64(ceil(fov_y / res_y)))

    class, loc = load_phantom_mat(objbrain; axis=axis, ss=ss, location=location)
    T1, T2, T2s, ρ = SpinProperty_1p5T(class)

    M, N = size(class)     # matrix size of the phantom mat

    if key == :ρ
        img = ρ
    elseif key == :T2
        img = T2
    elseif key == :T2s
        img = T2s
    elseif key == :T1
        img = T1
    elseif key == :raw
        img = class
    elseif key == :headmask
        img = (class.>0)*1 #All
    elseif key == :brainmask
        img = 
        (class.==23)*1 .+ #CSF
        (class.==46)*1 .+ #GM
        (class.==70)*1 .+ #WM
        (class.==93)*0 .+ #FAT1
        (class.==116)*0 .+ #MUSCLE
        (class.==139)*0 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*0 .+ #FAT2
        (class.==232)*0 .+ #DURA
        (class.==255)*0 #MARROW
    elseif key == :Δw
        if db0_type == :real
            B0map = load_B0map(db0_file; axis=axis, ss=1, location=location)
            img = imresize(B0map, size(class))
        elseif db0_type == :fat
            Δw_fat = γ * 1.5 * (-3.45) * 1e-6  # Hz
            img = (class.==93)*Δw_fat .+ #FAT1 
                (class.==209)*Δw_fat    #FAT2
        elseif db0_type == :quadratic
            img = quadraticFieldmap(size(class)...,db0_max)[:,:,1]
        end
    elseif key == :csm
        img = load_csm(csm_type, M, N, csm_nCoil; overlap=csm_overlap, 
            relative_radius=csm_radius, nRow=csm_nRow, nCol=csm_nCol, nBlock=csm_nBlock, verbose=verbose)
    end


    @info "PhantomReference" key=key axis=axis location=location obj_size=(M, N) center_range=center_range target_fov=target_fov target_size=target_size
    # img = img[get_center_range(M, center_range[1]), get_center_range(N, center_range[2])]
    img = get_center_crop(img, center_range[1], center_range[2])
    img = imresize(img, target_size)
    return Matrix(img')
end
