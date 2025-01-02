using MAT
import Statistics: quantile
using PyPlot, PyCall
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
mcolors = matplotlib.colors

mat_file = "$(@__DIR__)/demo/MultiChannel/snrInf_cgnr_40_L2_0.0.mat"   # *.mat file path
results = MAT.matread(mat_file);
#=
"csm"      => [a.u.] Coil-Sensitivity
"b0map"    => [Hz] B0 map
"headmask" => Head mask
"x_ref"    => Proton density map
"imgs"     => [a.u.] Image data, reconstruction results
"signal"   => [a.u.] Signal data, simulated signal
"labels"   => Image labels
=#
csm      = results["csm"];
B0map    = results["b0map"];
headmask = results["headmask"];
x_ref    = results["x_ref"];
imgs     = results["imgs"];
signal   = results["signal"];
labels   = results["labels"];

x, y = size(headmask);
headcountour = zeros(x, y);
α = zeros(x, y);
for i = 1:x, j = 1:x
    if headmask[i,j] == 1 
        if i == 1 || i == x || j == 1 || j == y
            headcountour[i,j] = 1
            α[i,j] = 1
        elseif (headmask[i-1,j] == 0 || headmask[i+1,j] == 0 || headmask[i,j-1] == 0 || headmask[i,j+1] == 0)
            headcountour[i,j] = 1
            α[i,j] = 1
        end
    end
end
img_headcountour = cat(headcountour, headcountour, headcountour, α; dims=3);  # with alpha channel
img_headref = cat(x_ref, x_ref, x_ref, headmask; dims=3);  # with alpha channel



matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=100)
# matplotlib.rc("font", family="Arial")    # "Times New Roman"
matplotlib.rcParams["mathtext.default"]
figure_width       = 5/2.54
figure_height      = 5/2.54
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

vmaxp              = 99;
vminp              = 1;



#############################################################################
# plot reconstructed images
#############################################################################
nImg, nX, nY = size(imgs)
# vmin, vmax = quantile(abs.(imgs)[:], 0), quantile(abs.(imgs)[:], 0.98)
for idx = 1 : nImg
    img = abs.(imgs)[idx,:,:]
    label = labels[idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    ax.set_title(label, fontsize=fontsize_label, color=color_label)
    fig.tight_layout(pad=0)
    # fig.savefig("$(outpath)/img_$(idx).png", dpi=300, transparent=true)
end

#############################################################################
# plot off-resonance map
#############################################################################
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1,1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end

ai = ax.imshow(B0map, cmap="jet")#, interpolation="gaussian") # "bilinear", "spline36", "gaussian"
# ax.imshow(img_headref, cmap="gray", alpha=0.1)
ax.imshow(img_headcountour, cmap="gray", alpha=1)
ax.set_title("ΔB₀ map", fontsize=fontsize_label, color=color_label)
divider = mpl_axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(ai, cax=cax)
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
cb.outline.set_visible(false)
cb.update_ticks()
fig.tight_layout(pad=0)

########################################################################
# plot proton density map
########################################################################
# colorbar setting for the ρ map of the phantom
ρ_values = sort(unique(x_ref))
cm_gray = plt.cm.get_cmap("gray")
cm_gray_ref = cm_gray.from_list("gray_ref",  cm_gray(ρ_values), length(ρ_values))
norm = mcolors.BoundaryNorm([-1; ρ_values] , cm_gray_ref.N)


fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1,1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end
ax.imshow(x_ref, cmap="gray") # "bilinear", "spline36", "gaussian"
divider = mpl_axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.02)
ax.set_title("Proton density", fontsize=fontsize_label, color=color_label)
cax.set_title("[a.u.]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm_gray_ref, norm=norm), cax=cax, drawedges=true)
cb.minorticks_off()
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel, length=ticklength, width=linewidth, pad=pad_labeltick)
cb.set_ticks(ρ_values - (ρ_values - [-1; ρ_values][1:end-1]) / 2)
cb.set_ticklabels([string(round(v, digits=2)) for v in ρ_values])
cb.ax.tick_params(size=0)
# cb.outline.set_visible(false)
cb.outline.set_edgecolor(color_label)
cb.outline.set_linewidth(linewidth)
cb.dividers.set_color(color_label)
cb.dividers.set_linewidth(linewidth)
fig.tight_layout(pad=0)

#############################################################################
# plot signal in each channel
#############################################################################
figure_width       = 12/2.54
figure_height      = 6/2.54

nSample, nCha = size(signal);
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_color(color_label)
    ax.spines[spine].set_visible(false)
end
ax.set_facecolor(color_facecolor)
fig.tight_layout(pad=0)
for cha = 1 : nCha
    ax.plot(abs.(signal[:, cha]), linewidth=0.5, color="C$(cha%9)", label="Channel $(cha)")
end
ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=3, 
    loc="upper right", bbox_to_anchor=(1,1),
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
fig.tight_layout(pad=0)
