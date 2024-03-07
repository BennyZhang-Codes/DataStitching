import KomaMRI.KomaMRIBase: brain_phantom2D
# import Base.show
abstract type PhantomType end
struct brain3D_02 <: PhantomType end

export brain3D_02


function info(s::Phantom)
	print("Phantom[name = $(s.name) | spins = $(length(s.x)) | x = $(minimum(s.x)*1e2):$(maximum(s.x)*1e2) cm | y = $(minimum(s.y)*1e2):$(maximum(s.y)*1e2) cm | z = $(minimum(s.z)*1e2):$(maximum(s.z)*1e2) cm ]")
end


"""
    obj = brain_phantom2D(p::brain3D_02, location::AbstractFloat; axis="axial", ss=4)

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

# Returns
- `obj`: (`::Phantom`) Phantom struct

# Examples
```julia-repl
julia> obj = brain_phantom2D(p::brain3D_02, location::AbstractFloat; axis="sagittal", ss=1)

julia> plot_phantom_map(obj, :ρ)
```
"""
function brain_phantom2D(p::brain3D_02; axis="axial", ss=4, location=0.5) :: Phantom
    @assert 0 <= location <= 1
    # Get data from .mat file
    path = @__DIR__
    data = MAT.matread(path*"/brain3D_0.2.mat")["data"]
    
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

    # Define spin position vectors
    Δx = .2e-3*ss
    M, N = size(class)
    FOVx = (M-1)*Δx #[m]
    FOVy = (N-1)*Δx #[m]
    x = -FOVx/2:Δx:FOVx/2 #spin coordinates
    y = -FOVy/2:Δx:FOVy/2 #spin coordinates
    x, y = x .+ y'*0, x*0 .+ y' #grid points

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
        name = "brain2D_02_"*axis,
		x = y[ρ.!=0],
		y = x[ρ.!=0],
		z = 0*x[ρ.!=0],
		ρ = ρ[ρ.!=0],
		T1 = T1[ρ.!=0],
		T2 = T2[ρ.!=0],
		T2s = T2s[ρ.!=0],
		Δw = Δw[ρ.!=0],
    )
	return obj
end
