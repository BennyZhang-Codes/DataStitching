using HighOrderMRI, PyPlot, MAT
import Statistics: quantile
using PyPlot, PyCall
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")

T = Float64


path     = "E:/pulseq/20250222_LYQ"
outpath = "$(@__DIR__)/Figures/FigS2/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

gre_file = "$(path)/gre/gre_H15_1p0_200_r1.mat"
@info "gre file: $(gre_file)"

# gre data
gre_data = matread(gre_file);
gre_img  = gre_data["img"];
gre_csm  = gre_data["csm"];
gre_b0   = gre_data["b0"];
gre_mask = gre_data["mask"];


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
# matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 5/2.53999863
figure_height      = 5/2.53999863
linewidth          = 0.5
ticklength         = 1.5

fontsize_subfigure = 9
fontsize_label     = 7
fontsize_legend    = 7
fontsize_ticklabel = 6

pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#ffffff"
color_label        = "#000000"

b0 = gre_b0
fig = plt.figure(figsize=(figure_width, figure_height), facecolor=color_facecolor)
ax    = fig.add_axes([0, 0, 0.7, 0.7])
cax   = fig.add_axes([0.72, 0, 0.05*0.7, 0.7])

ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end

ai = ax.imshow(b0, cmap="jet", vmin=-120, vmax=120)#, interpolation="gaussian") # "bilinear", "spline36", "gaussian"
# divider = mpl_axes_grid1.make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(ai, cax=cax)
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
cb.outline.set_visible(false)
# cb.set_ticks([-120.0, -60.0, 0.0, 60.0, 120.0])
cb.update_ticks()
fig.tight_layout(pad=0)


fig.savefig("$(outpath)/B0map.png", dpi=900, transparent=true, bbox_inches="tight", pad_inches=0)
fig.savefig("$(outpath)/B0map.svg", dpi=900, transparent=true, bbox_inches="tight", pad_inches=0.)

