using PyPlot
using MAT
import Statistics: quantile

dir = "Figures/out/recon_fully_sampled"; if ispath(dir) == false mkpath(dir) end     # output directory

matfile = "recon_fully_sampled_admm_35_L1_1e-3"

imgs = MAT.matread("$(dir)/$(matfile).mat")["imgs"];
nFrame, nX, nY = size(imgs)

figure_width       = 10
figure_height      = 3
vmaxp              = 99
vminp              = 1 
cmap               = "gray"
fontsize_legend    = 10
fontsize_label     = 12
fontsize_ticklabel = 8
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_imagelabel   = "#ffffff"

fig, axs = plt.subplots(nrows=1, ncols=nFrame, figsize=(figure_width, figure_height), facecolor=color_facecoler)

for idx = 1 : nFrame
    ax = axs[idx]
    img = imgs[idx,:,:]
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    ax.axis("off")
    ax.text(0.02, 0.98, "$(Char(96+idx))", fontsize=fontsize_label, color=color_imagelabel, transform=ax.transAxes, ha="left", va="top")
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
end
axs[1].set_title("w/o ΔB₀\nnominal"  , fontsize=fontsize_label, color=color_label)
axs[2].set_title("w/ ΔB₀\nnominal"      , fontsize=fontsize_label, color=color_label)
axs[3].set_title("w/ ΔB₀\nstitching 110", fontsize=fontsize_label, color=color_label)
axs[4].set_title("w/ ΔB₀\nstitching 111", fontsize=fontsize_label, color=color_label)
axs[5].set_title("w/ ΔB₀\nconventional 111" , fontsize=fontsize_label, color=color_label)

fig.tight_layout(pad=0, w_pad=0.3, h_pad=0)

fig.savefig("$(dir)/$(matfile).png", dpi=300, bbox_inches="tight")





############

imgs_cgnr = MAT.matread("$(dir)/recon_fully_sampled_cgnr_100_L2_1e-9.mat")["imgs"];
imgs_admm = MAT.matread("$(dir)/recon_fully_sampled_admm_20_TV_1e-4.mat")["imgs"];
imgs = cat(imgs_cgnr, imgs_admm, dims=1)
nFrame, nX, nY = size(imgs)
nRow = 2
nCol = nFrame ÷ nRow

figure_width       = 10
figure_height      = 5
vmaxp              = 99
vminp              = 1      
cmap               = "gray"
fontsize_legend    = 10
fontsize_label     = 14
fontsize_imagelabel= 13
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_imagelabel   = "#ffffff"

fig, axs = plt.subplots(nrows=nRow, ncols=nCol, figsize=(figure_width, figure_height), facecolor=color_facecoler)

for idx = 1 : nFrame
    row, col = Int(ceil(idx/nCol)), (idx-1)%nCol + 1
    ax = axs[row, col]
    img = imgs[idx,:,:]
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.tick_params(axis="both", color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.text(0.02, 0.98, "$(Char(96+idx))", fontsize=fontsize_imagelabel, color=color_imagelabel, transform=ax.transAxes, ha="left", va="top")
end
axs[1,1].set_ylabel("CG"  , fontsize=fontsize_label, color=color_label)
axs[2,1].set_ylabel("ADMM", fontsize=fontsize_label, color=color_label)
axs[1,1].set_title("w/o ΔB₀\nnominal"  , fontsize=fontsize_label, color=color_label)
axs[1,2].set_title("w/ ΔB₀\nnominal"      , fontsize=fontsize_label, color=color_label)
axs[1,3].set_title("w/ ΔB₀\nstitching 110", fontsize=fontsize_label, color=color_label)
axs[1,4].set_title("w/ ΔB₀\nstitching 111", fontsize=fontsize_label, color=color_label)
axs[1,5].set_title("w/ ΔB₀\nconventional 111" , fontsize=fontsize_label, color=color_label)

fig.tight_layout(w_pad=0, h_pad=0, pad=0)

fig.savefig("$(dir)/recon_fully_sampled.png", dpi=300, bbox_inches="tight")