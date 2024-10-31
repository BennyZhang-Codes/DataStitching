include("Figures/Fig_preset.jl");
using MAT
import Statistics: quantile

dir = "workplace/Abstract/Simulation/out"; if ispath(dir) == false mkpath(dir) end     # output directory

###### load images reconstructed by different considerations of the extended signal model 
solver = "admm";
regularization = "TV";
λ = 1.e-4;
iter=20;
matfile = "fully_$(solver)_$(iter)_$(regularization)_$(λ)"
imgs = MAT.matread("$(dir)/$(matfile).mat")["imgs"];
nFrame, nX, nY = size(imgs)

###### get the reference images
simtype = SimType(B0=true, T2=false, ss=5);
phantom = BrainPhantom(prefix="brain3D", x=0.2, y=0.2, z=0.2); # decide which phantom file to use
maxOffresonance = 200.;        



figure_width       = 5/2.54
figure_height      = 5/2.54
vmaxp              = 99;
vminp              = 1;


########################################################################
# Coil-Sensitivity Map
########################################################################
csm = KomaHighOrder.csm_Birdcage(215, 180, 8) 

vmin, vmax = quantile(abs.(csm)[:], 0), quantile(abs.(csm)[:], 1)
for idx = 1 : 8
    img = abs.(csm)[:,:, idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("workplace/Abstract/Figure1/out/Sim_csm_Birdcage_cha$(idx).png", dpi=300, transparent=true, bbox_inches="tight", pad_inches=0)
end



########################################################################
# Off-resonance map
########################################################################
# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=1, location=0.5, target_fov=(215, 180), target_resolution=(1,1), B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
B0map = rotl90(B0map);

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1,1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end

ai = ax.imshow(B0map, cmap="jet")#, interpolation="gaussian") # "bilinear", "spline36", "gaussian"
fig.tight_layout(pad=0)
fig.savefig("workplace/Abstract/Figure1/out/Sim_b0_quadratic.png", dpi=300, transparent=true, bbox_inches="tight", pad_inches=0)


########################################################################
# Proton density map
########################################################################
x_ref = brain_phantom2D_reference(phantom; ss=1, location=0.5, key=:ρ, target_fov=(215, 180), target_resolution=(1,1));
x_ref = rotl90(x_ref);

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1,1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end
ax.imshow(x_ref, cmap="gray") # "bilinear", "spline36", "gaussian"
fig.tight_layout(pad=0)
fig.savefig("workplace/Abstract/Figure1/out/Sim_ProtonDensity.png", dpi=300, transparent=true, bbox_inches="tight", pad_inches=0)





