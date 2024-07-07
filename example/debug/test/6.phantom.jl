location = 0.5
ss = 4
data = MAT.matread("src\\phantom\\brain3D_0.2.mat")["data"]

M, N, Z = size(data)

z = Int32(ceil(N*location))
class = data[1:ss:end,z,1:ss:end]

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
    name = "brain2D_02",
    x = y[ρ.!=0],
    y = x[ρ.!=0],
    z = 0*x[ρ.!=0],
    ρ = ρ[ρ.!=0],
    T1 = T1[ρ.!=0],
    T2 = T2[ρ.!=0],
    T2s = T2s[ρ.!=0],
    Δw = Δw[ρ.!=0],
)

plot_phantom_map(obj, :ρ)


# h = scatter3d(x=obj.x, y=obj.y, z=obj.z)
# plot(h)


using MAT
path = "src/phantom"
data = MAT.matread(path*"/brain3D_(400, 362, 434)_(0.025, 0.5, 0.5).mat")["data"]
ss = 9

class = data[1:ss:end,1:ss:end,1:1:400]
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
    name = "brain3D_Z0.025",
    x = y[ρ.!=0],
    y = x[ρ.!=0],
    z = z[ρ.!=0],
    ρ = ρ[ρ.!=0],
    T1 = T1[ρ.!=0],
    T2 = T2[ρ.!=0],
    T2s = T2s[ρ.!=0],
    Δw = Δw[ρ.!=0],
)
info(obj)
plot_phantom_map(obj, :ρ)

