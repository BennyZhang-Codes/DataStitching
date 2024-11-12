include("Figures/Fig_preset.jl");
using MAT
import Statistics: quantile

dir = "workplace/Abstract/Simulation/out"; if ispath(dir) == false mkpath(dir) end     # output directory
outpath = "workplace/Abstract/Figures/Figure2/out"; if ispath(outpath) == false mkpath(outpath) end


###### get the reference images
Nx = Ny = 150;

csmtype= :birdcage; nCoil=8; nrows=2; ncols=4;
coil = csm_Birdcage(217, 181, nCoil, verbose=true);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
csm = permutedims(mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]), [3, 1, 2]);
fig_csm = plt_images(csm,width=8, height=8)

figure_width       = 5/2.54
figure_height      = 5/2.54
vmaxp              = 99;
vminp              = 1;


########################################################################
# Coil-Sensitivity Map
########################################################################
vmin, vmax = quantile(abs.(csm)[:], 0), quantile(abs.(csm)[:], 1)
for idx = 1 : 8
    img = abs.(csm)[idx,:,:]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/Sim_csm_Birdcage_cha$(idx).png", dpi=900, transparent=true, bbox_inches="tight", pad_inches=0)
end




