

function brain_phantom2D_reference(
    objbrain::BrainPhantom; 
    axis="axial", 
    ss::Int64=3, 
    location::Float64=0.8, 
    key::Symbol=:ρ, 
    B0_type::Symbol=:fat,          # load B0 map
    B0_file::Symbol=:B0,
    maxOffresonance::Float64=125., # max off-resonance
    target_fov=(150, 150), 
    target_resolution=(1,1))

    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert key in [:ρ, :T2, :T2s, :T1, :Δw, :raw, :headmask, :brainmask] "key must be ρ, T2, T2s, T1, Δw, raw, headmask or brainmask"
    @assert B0_type in [:real, :fat, :quadratic] "B0_type must be one of the following: :real, :fat, :quadratic"
    @assert 0 <= location <= 1 "location must be between 0 and 1"

    Δx = Δy = 0.2   # resolution of phantom: phantom_dict[:brain2d]
    fov_x, fov_y = target_fov
    res_x, res_y = target_resolution
    center_range = (Int64(ceil(fov_x / (Δx * ss))), Int64(ceil(fov_y / (Δy * ss)))) 
    target_size = (Int64(ceil(fov_x / res_x)), Int64(ceil(fov_y / res_y)))

    data = MAT.matread(objbrain.matpath)["data"]

    M, N, Z = size(data)
    if axis == "axial"
        z = Int32(ceil(Z*location))
        class = data[1:ss:end,1:ss:end, z]
    elseif axis == "coronal"
        m = Int32(ceil(M*location))
        class = data[m, 1:ss:end,1:ss:end]   
    elseif axis == "sagittal"
        n = Int32(ceil(N*location))
        class = data[1:ss:end, n,1:ss:end]
    end
    if key == :ρ
        img = (class.==23)*1 .+ #CSF
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
    elseif key == :T2
        img = (class.==23)*329 .+ #CSF
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
        img = img .* 1e-3  # s
    elseif key == :T2s
        img = (class.==23)*58 .+ #CSF
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
        img = img .* 1e-3  # s
    elseif key == :T1
        img = (class.==23)*2569 .+ #CSF
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
        img = img .* 1e-3  # s
    elseif key == :Δw
        if B0_type == :real
            B0map = brain_phantom2D_B0map(B0_file; axis=axis, ss=1, location=location)
            img = imresize(B0map, size(class))
        elseif B0_type == :fat
            Δw_fat = γ * 1.5 * (-3.45) * 1e-6  # Hz
            img = (class.==93)*Δw_fat .+ #FAT1 
                (class.==209)*Δw_fat    #FAT2
        elseif B0_type == :quadratic
            img = quadraticFieldmap(size(class)...,maxOffresonance)[:,:,1]
        end
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
    end
    M, N = size(img)

    @info "PhantomReference" key=key axis=axis location=location obj_size=(M, N) center_range=center_range target_fov=target_fov target_size=target_size
    img = img[get_center_range(M, center_range[1]), get_center_range(N, center_range[2])]
    img = imresize(img, target_size)
    return Matrix(img')
end
