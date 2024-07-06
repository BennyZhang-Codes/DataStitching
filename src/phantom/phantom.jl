# definition of the abstract PhantomType and BrainPhantom struct
include("PhantomType.jl")
export BrainPhantom

include("hohantom2d.jl")

# include the Coil-Sensitivity Map (CSM) module
include("csm/csm.jl")

include("B0map/B0map.jl")



# functions to generate different types of phantoms
include("phantom2d.jl")
include("phantom2d_reference.jl")
include("phantom3d.jl")

export brain_phantom2D_reference

# function to print the information of a Phantom object
function info(s::Phantom)
	print("Phantom[name = $(s.name) | spins = $(length(s.x)) | x = $(minimum(s.x)*1e2):$(maximum(s.x)*1e2) cm | y = $(minimum(s.y)*1e2):$(maximum(s.y)*1e2) cm | z = $(minimum(s.z)*1e2):$(maximum(s.z)*1e2) cm ]")
    print("\n")
end
export info



##########
# obj = brain_phantom2D(brain2D(); ss=3, location=0.8, B0map=:quadratic, maxOffresonance=10.); info(obj); plot_phantom_map(obj, :Δw)
# ref = brain_phantom2D_reference(brain2D();B0map=:quadratic,key=:Δw, maxOffresonance=10.); plot_image(ref,zmin=-10)

# B0map = brain_phantom2D_reference(brain2D(); ss=3, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
#                                            B0map=:quadratic,key=:Δw, maxOffresonance=5.);
# Nx = Ny = 150
# Δx = Δy = 1e-3 # m

# x, y = 1:Nx, 1:Ny
# x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
# x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1

# plot_bgcolor = "rgb(22,26,29)"
# grid_color = "rgb(40,52,66)"
# p1 = plot(scatter3d(x=x, y=y, z=vec(B0map),       marker=attr(size=0.7), mode="markers"), Layout(
#     paper_bgcolor="rgba(0,0,0,0)", 
#     scene=attr(xaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
#                 yaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
#                 zaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color)),
#     font=attr(family="Times New Roman",color="gray"),
#     width=500,height=500))

# savefig(p1, "/B0map3D.svg",            width=500,height=500,format="svg")