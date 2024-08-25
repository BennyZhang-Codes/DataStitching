using PyPlot
using MAT
import Statistics: quantile

dir = "Figures/Fig5"; if ispath(dir) == false mkpath(dir) end     # output directory
solver = "admm";
regularization = "TV";
λ = 1.e-4;
iter=20;

matfile = "fully_$(solver)_$(iter)_$(regularization)_$(λ)"

imgs = MAT.matread("$(dir)/$(matfile).mat")["imgs"];
nFrame, nX, nY = size(imgs)


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 12/2.54
figure_height      = 6/2.54
vmaxp              = 99
vminp              = 1 
cmap               = "gray"
fontsize_legend    = 5
fontsize_label     = 8
fontsize_subfigure = 8
pad_label          = 2
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_subfigure    = "#ffffff"

fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(figure_width, figure_height), facecolor=color_facecoler)

idx_img = [1 2 3 0; 4 5 6 7]
for row = 1 : 2
    for col = 1 : 4
        ax = axs[row, col]
        idx = idx_img[row, col]

        ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
        for spine in ax.spines  # "left", "right", "bottom", "top"
            ax.spines[spine].set_visible(false)
        end
        if idx == 0
            
            continue
        end
        img = imgs[idx,:,:]
        vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
        ax.text(0.02, 0.98, "$(Char(96+idx))", fontsize=fontsize_subfigure, color=color_subfigure, transform=ax.transAxes, ha="left", va="top")
        ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    end
end

axs[1,1].set_ylabel("w/o ΔB₀"         , fontsize=fontsize_label, color=color_label, labelpad=pad_label, rotation=0, ha="right", va="center", x=0, y=0.5,)
axs[2,1].set_ylabel("with ΔB₀"        , fontsize=fontsize_label, color=color_label, labelpad=pad_label, rotation=0, ha="right", va="center", x=0, y=0.5,)

axs[1,1].set_title("nominal"          , fontsize=fontsize_label, color=color_label)
axs[1,2].set_title("stitching 110"    , fontsize=fontsize_label, color=color_label)
axs[1,3].set_title("stitching 111"    , fontsize=fontsize_label, color=color_label)
axs[2,4].set_title("conventional 111" , fontsize=fontsize_label, color=color_label)

# fig.tight_layout(pad=0, w_pad=0, h_pad=0)
fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0, hspace=0)
fig.savefig("$(dir)/Fig5_$(matfile).png", dpi=300, bbox_inches="tight")





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