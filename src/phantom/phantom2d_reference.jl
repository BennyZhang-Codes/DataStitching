
function get_center_range(x::Int64, x_range::Int64)
    center = Int64(floor(x/2))
    return center - Int64(ceil(x_range/2))+1 : center + 1 + Int64(ceil(x_range/2))-1
end


function brain_phantom2D_reference(p::PhantomType; axis="axial", ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1))
    path = (@__DIR__) * phantom_dict[:path]
    @assert axis in ["axial", "coronal", "sagittal"] "axis must be one of the following: axial, coronal, sagittal"
    @assert key in [:ρ, :T2, :T2s, :T1, :Δw, :raw, :binary] "key must be ρ, T2, T2s, T1, Δw, raw or binary"
    @assert 0 <= location <= 1 "location must be between 0 and 1"

    @assert isfile(path*"/$(p.file)") "the phantom file does not exist: $(path*"/$(p.file)")"
    Δx = Δy = 0.2   # resolution of phantom: phantom_dict[:brain2d]
    fov_x, fov_y = target_fov
    res_x, res_y = target_resolution

    data = MAT.matread(path*"/$(p.file)")["data"]

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
    elseif key == :Δw
        B0map = brain_phantom2D_B0map(; axis=axis, ss=1, location=location)
        img = imresize(B0map, size(class))
    elseif key == :raw
        img = class
    elseif key == :binary
        img = (class.>0)*1 #All
    end
    M, N = size(img)
    center_range = (Int64(ceil(fov_x / (Δx * ss))), Int64(ceil(fov_y / (Δy * ss)))) 
    target_size = (Int64(ceil(fov_x / res_x)), Int64(ceil(fov_y / res_y)))
    @info "PhantomReference" key=key axis=axis location=location obj_size=(M, N) center_range=center_range target_fov=target_fov target_size=target_size
    img = img[get_center_range(M, center_range[1]), get_center_range(N, center_range[2])]
    img = imresize(img, target_size)
    return Matrix(img')
end
