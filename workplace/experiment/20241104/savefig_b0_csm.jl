using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData
import Statistics: quantile
using PyPlot, PyCall
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

T = Float64

# Coil-Sensitivity Map (CSM), ΔB0 map, mask
# b0   = -matread(gre_file)["b0"];
# csm  = matread(gre_file)["csm"];    # (nY, nX, nCha)
# mask = matread(gre_file)["mask"];
# nY, nX, nCha = size(csm);




path         = "E:/pulseq/20241104_ABDL/"
gre_seq_file = "$(path)/seq/gres6e_fov200_200_tr25_fa10_bw556.seq" 
gre_mrd_file = "$(path)/mrd/meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mrd"

outpath = "$(path)/out/gre"; if ispath(outpath) == false mkpath(outpath) end

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
fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=12/2.54, height=6/2.54)  # (nCha, nY, nX)
fig.savefig("$(outpath)/csm.png", dpi=900, transparent=true)
fig.savefig("$(outpath)/csm.svg", dpi=900, transparent=true)


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
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/csm_$(idx).png", dpi=900, transparent=true)
    fig.savefig("$(outpath)/csm_$(idx).svg", dpi=900, transparent=true)
end


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 6/2.54
figure_height      = 6/2.54
linewidth          = 0.5
ticklength         = 1.5
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 7
fontsize_subfigure = 8
pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#ffffff"
color_label        = "#000000"



############
# ΔB₀
############
TE = TE;  # s
smap = csm[:,:,:]; # [1,2,3,5,6,7,8] (nY, nX, nCha)
# create mask from Coil-Sensitivity Map
mask = get_mask(smap[:,:,1], threshold=0); #plt_image(mask)
images = convert_ifft(kdata[:,:,:,:], dims=[2,3]); # (nCha, nY, nX, nCon)
ydata = permutedims(images, [2,3,1,4]);            # (nCha, nY, nX, nCon) => (nY, nX, nCha, nCon)

b0, yik_sos = get_B0map(ydata, TE, smap, mask);


vmin, vmax = quantile(abs.(yik_sos)[:], 0), quantile(abs.(yik_sos)[:], 0.995)
for idx = 1 : nCon
    img = abs.(yik_sos)[:,:, idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/gre_mag_echo$(idx).png", dpi=900, transparent=true)
    fig.savefig("$(outpath)/gre_mag_echo$(idx).svg", dpi=900, transparent=true)
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
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/gre_pha_echo$(idx).png", dpi=900, transparent=true)
    fig.savefig("$(outpath)/gre_pha_echo$(idx).svg", dpi=900, transparent=true)
end

b0 = -b0


fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1,1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end

ai = ax.imshow(b0, cmap="jet", vmin=-150, vmax=150)#, interpolation="gaussian") # "bilinear", "spline36", "gaussian"
# ax.imshow(img_headref, cmap="gray", alpha=0.1)
# ax.imshow(img_headcountour, cmap="gray", alpha=1)
# ax.set_title("Off-resonance", fontsize=fontsize_label, color=color_label)
divider = mpl_axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(ai, cax=cax)
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
cb.outline.set_visible(false)
cb.update_ticks()
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/B0map.png", dpi=900, transparent=true)
fig.savefig("$(outpath)/B0map.svg", dpi=900, transparent=true)













####### 

function plt_images1(
    imgs::AbstractArray{<:Real, 3};
    nRow               = nothing  ,
    nCol               = nothing  ,
    title              = ""       ,
    width              = 5        ,
    height             = 5        ,
    vmaxp              = 100      ,
    vminp              = 0        ,
    cmap               = "gray"   ,
    fontsize_title     = 10       ,
    color_facecolor    = "#ffffff",
    )
    nFrame, nX, nY = size(imgs)
    if nRow === nothing || nCol === nothing
        nCol, nRow = get_factors(nFrame)
    end

    vmin, vmax = quantile(reshape(imgs, :), vminp/100), quantile(reshape(imgs, :), vmaxp/100)

    fig, axs = plt.subplots(nrows=nRow, ncols=nCol, figsize=(width, height), facecolor=color_facecolor)
	axs = reshape(axs, nRow, nCol) # to keep axs as a 2D array
    
    fig.suptitle(title, fontsize=fontsize_title)
    for frame = 1:nFrame
        row, col = Int(ceil(frame/nCol)), (frame-1)%nCol + 1
        ax = axs[row, col]
        ax.imshow(imgs[frame, :, :], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.axis("off")
    end
    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    return fig
end

fig = plt_images1(permutedims(abs.(csm), [3,1,2]),width=6/2.54, height=12/2.54)  # (nCha, nY, nX)
fig.savefig("$(outpath)/csm_vertical.png", dpi=900, transparent=true)
fig.savefig("$(outpath)/csm_vertical.svg", dpi=900, transparent=true)