import KomaMRI.KomaMRIBase: brain_phantom3D
# import Base.show

Base.@kwdef mutable struct brain3D <: PhantomType 
    name::String = "3D phantom with 0.025 mm super-resolution in Z-axis"
    file::String = "(171, 191)_brain3D_(400, 362, 434)_(0.025, 0.5, 0.5).mat"
end


export brain3D


"""
    obj = brain_phantom3D(p::brain3D; ss=4, start_end=[150, 250])

Creates a three-dimensional brain Phantom struct.

# Returns
- `obj`: (`::Phantom`) Phantom struct
```
"""
function brain_phantom3D(p::PhantomType; ss=4, start_end=[150, 250]) :: Phantom
    
    # Get data from .mat file
    path = @__DIR__
    @assert isfile(path*"/$(p.file)") "File $(p.file) does not exist in $(path)"
    data = MAT.matread(path*"/$(p.file)")["data"]
    _,_,Nz = size(data)
    @assert start_end[1] >= 1 && start_end[2] <= Nz "start_end must be between 1 and $Nz"


    class = data[1:ss:end,1:ss:end,start_end[1]:1:start_end[2]]

    # Define spin position vectors
    Δx = Δy = .5e-3*ss
    Δz= .025e-3
    M, N, Z = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δy #[m]
    FOVz = (Z-1)*Δz #[m]
    xx = reshape(-FOVx/2:Δx:FOVx/2,M,1,1) #spin coordinates
    yy = reshape(-FOVy/2:Δy:FOVy/2,1,N,1) #spin coordinates
    zz = reshape(-FOVz/2:Δz:FOVz/2,1,1,Z) #spin coordinates
    x = 1*xx .+ 0*yy .+ 0*zz
    y = 0*xx .+ 1*yy .+ 0*zz
    z = 0*xx .+ 0*yy .+ 1*zz

    # Define spin property vectors
    T2 = (class.==23)*329 .+ #CSF
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
    T2s = (class.==23)*58 .+ #CSF
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
    T1 = (class.==23)*2569 .+ #CSF
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
    ρ = (class.==23)*1 .+ #CSF
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
	Δw_fat = -220*2π
	Δw = (class.==93)*Δw_fat .+ #FAT1
		(class.==209)*Δw_fat    #FAT2
	T1 = T1*1e-3
	T2 = T2*1e-3
	T2s = T2s*1e-3

    # Define and return the Phantom struct
    obj = Phantom{Float64}(
        name = "brain3D_ss$(ss)_[$(start_end[1]), $(start_end[2])]",
		x = y[ρ.!=0],
		y = x[ρ.!=0],
		z = z[ρ.!=0],
		ρ = ρ[ρ.!=0],
		T1 = T1[ρ.!=0],
		T2 = T2[ρ.!=0],
		T2s = T2s[ρ.!=0],
		Δw = Δw[ρ.!=0],
    )
	return obj
end
