using PyPlot
using MAT
import Statistics: quantile

dir = "Figures/out"; if ispath(dir) == false mkpath(dir) end     # output directory

matfile = "recon_fully_sampled_admm_50_L1_1e-3"

imgs = MAT.matread("$(dir)/$(matfile).mat")["imgs"];
nFrame, nX, nY = size(imgs)

figure_width       = 10
figure_height      = 2
vmaxp              = 100 
vminp              = 0       
cmap               = "gray"
fontsize_legend    = 10
fontsize_label     = 12
fontsize_ticklabel = 8
color_facecoler    = "#ffffff"
color_label        = "#ffffff"

fig, axs = plt.subplots(nrows=1, ncols=nFrame, figsize=(figure_width, figure_height), facecolor=color_facecoler)

for idx = 1 : nFrame
    ax = axs[idx]
    img = imgs[idx,:,:]
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    ax.axis("off")
    ax.text(0.02, 0.98, "$(Char(96+idx))", fontsize=fontsize_label, color=color_label, transform=ax.transAxes, ha="left", va="top")
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
end
fig.tight_layout(pad=0, w_pad=0.3, h_pad=0)

fig.savefig("$(dir)/$(matfile).png", dpi=300, bbox_inches="tight")