

function brain_phantom2D_reference(
    objbrain::BrainPhantom; 
    axis="axial", 
    ss::Int64=3, 
    location::Float64=0.8, 
    key::Symbol=:ρ, 
    B0type::Symbol=:fat,           # load B0 map
    B0_file::Symbol=:B0,
    maxOffresonance::Float64=125., # max off-resonance
    target_fov=(150, 150), 
    target_resolution=(1,1))

    @assert key in [:ρ, :T2, :T2s, :T1, :Δw, :raw, :headmask, :brainmask] "key must be ρ, T2, T2s, T1, Δw, raw, headmask or brainmask"
    @assert B0type in [:real, :fat, :quadratic] "B0_type must be one of the following: :real, :fat, :quadratic"
    @assert 0 <= location <= 1 "location must be between 0 and 1"

    Δx = Δy = 0.2   # resolution of phantom: phantom_dict[:brain2d]
    fov_x, fov_y = target_fov
    res_x, res_y = target_resolution
    center_range = (Int64(ceil(fov_x / (Δx * ss))), Int64(ceil(fov_y / (Δy * ss)))) 
    target_size = (Int64(ceil(fov_x / res_x)), Int64(ceil(fov_y / res_y)))

    class = load_phantom_mat(objbrain; axis=axis, ss=ss, location=location)
    T1, T2, T2s, ρ = SpinProperty_1p5T(class)

    if key == :ρ
        img = ρ
    elseif key == :T2
        img = T2
    elseif key == :T2s
        img = T2s
    elseif key == :T1
        img = T1
    elseif key == :Δw
        if B0type == :real
            B0map = load_B0map(B0_file; axis=axis, ss=1, location=location)
            img = imresize(B0map, size(class))
        elseif B0type == :fat
            Δw_fat = γ * 1.5 * (-3.45) * 1e-6  # Hz
            img = (class.==93)*Δw_fat .+ #FAT1 
                (class.==209)*Δw_fat    #FAT2
        elseif B0type == :quadratic
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
