using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

#=
prepare data for synchronization
1. Multi-Echo GRE (ME-GRE) data: 
    (1) Coil-Sensitivity Map (CSM)
    (2) ΔB₀ map
2. Spiral imaging datas 
    (1) MRI Signal data

Notes: directory structure:
    - path/seq/*.seq     # sequence files
    - path/mrd/*.mrd     # MRD files
    - path/syn/*.mat     # output directory
=#
figure_width       = 5/2.54
figure_height      = 5/2.54
vmaxp              = 99;
vminp              = 1;


path         = "E:/pulseq/20241010_skope_fa90/invivo"
gre_seq_file = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("gres6", f)][1]
gre_mrd_file = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin("gres6", f)][1]

outpath = "workplace/Abstract/Figures/Figure1/out"; if ispath(outpath) == false mkpath(outpath) end

kdata, ktraj, kdims, shape, TE, ReadoutMode, acqData, raw = read_gre(gre_seq_file, gre_mrd_file);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;

# IFFT of gre
gre_imgs = convert_ifft(kdata, dims=[2,3]);
gre_imgs = CoilCombineSOS(abs.(gre_imgs), 1);
gre_imgs = permutedims(gre_imgs, [3,1,2]);
fig = plt_images(gre_imgs, width=7.5, height=5, vminp=0, vmaxp=99)


############
# espirit
############
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)
csm = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nX, nY, 1, nCha) => (nY, nX, nCha, 1) => (nY, nX, nCha)
fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=10, height=5)  # (nCha, nY, nX)
fig.savefig("$(path)/syn/$(basename(gre_mrd_file)[1:end-4])_CoilSens_espirit.png", dpi=300, bbox_inches="tight", pad_inches=0)

vmin, vmax = quantile(abs.(csm)[:], 0), quantile(abs.(csm)[:], 1)
for idx = 1 : nCha
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
    fig.savefig("$(outpath)/csm_$(idx).png", dpi=300, transparent=true)
end


############
# ΔB₀
############
TE = TE;  # s
smap = csm[:,:,:]; # [1,2,3,5,6,7,8] (nY, nX, nCha)
# create mask from Coil-Sensitivity Map
mask = get_mask(smap[:,:,1], threshold=0); plt_image(mask)
images = convert_ifft(kdata[:,:,:,:], dims=[2,3]); # (nCha, nY, nX, nCon)
ydata = permutedims(images, [2,3,1,4]);            # (nCha, nY, nX, nCon) => (nY, nX, nCha, nCon)

b0, yik_sos = get_B0map(ydata, TE, smap, mask);
fig = plt_images(permutedims(  abs.(yik_sos), [3,1,2]),width=7.5, height=5)
fig = plt_images(permutedims(angle.(yik_sos), [3,1,2]),width=7.5, height=5)
fig = plt_B0map(b0, width=7.5, height=5)

vmin, vmax = quantile(abs.(yik_sos)[:], 0), quantile(abs.(yik_sos)[:], 0.99)
for idx = 1 : nCon
    img = abs.(yik_sos)[:,:, idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/gre_mag_echo$(idx).png", dpi=300, transparent=true)
end

vmin, vmax = quantile(angle.(yik_sos)[:], 0), quantile(angle.(yik_sos)[:], 1)
for idx = 1 : nCon
    img = angle.(yik_sos)[:,:, idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/gre_pha_echo$(idx).png", dpi=300, transparent=true)
end

img = -b0

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1, 1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end
ax.imshow(img, cmap="jet", vmin=-130,   vmax=130)
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/B0map.png", dpi=300, transparent=false)