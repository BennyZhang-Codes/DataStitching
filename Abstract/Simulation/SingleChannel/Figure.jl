include("Figures/Fig_preset.jl");
using MAT
import Statistics: quantile

outpath = "$(@__DIR__)/Abstract/Simulation/SingleChannel/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory
###### load images reconstructed by different considerations of the extended signal model 
solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;
matfile = "fully_$(solver)_$(iter)_$(regularization)_$(λ)"
mat = MAT.matread("$(outpath)/$(matfile).mat");

"imgs"=>imgs, "lables"=>lables, "B0map"=>B0map, "x_ref"=>x_ref, "headmask"=>headmask
imgs   = mat["imgs"]
labels = mat["labels"]
B0map  = mat["B0map"]
x_ref  = mat["x_ref"]
headmask = mat["headmask"]
nFrame, nX, nY = size(imgs)

###### get the reference images
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


using PyPlot, PyCall
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
mtransforms = PyPlot.plt.matplotlib.transforms
mcolors = matplotlib.colors

matplotlib.rc("mathtext", default="regular");
matplotlib.rc("figure", dpi=200);
matplotlib.rc("font", family="Arial");
matplotlib.rcParams["mathtext.default"];
figure_width       = 5/2.54;
figure_height      = 5/2.54;
vmaxp              = 99;
vminp              = 1;
cmap               = "gray";
fontsize_legend    = 5;
fontsize_label     = 6;
fontsize_subfigure = 8;
fontsize_ticklabel = 4;
linewidth          = 0.5;
ticklength         = 1.5;
pad_label          = 2;
pad_labeltick      = 2;
color_facecolor    = "#ffffff";
color_label        = "#000000";
color_subfigure    = "#ffffff";

vmaxp              = 99;
vminp              = 1;

########################################################################
# Off-resonance map
########################################################################
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
ax.set_title("Off-resonance", fontsize=fontsize_label, color=color_label)
divider = mpl_axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(ai, cax=cax)
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
cb.outline.set_visible(false)
cb.update_ticks()
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/B0map.png", dpi=300, transparent=true)


########################################################################
# Proton density map
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
fig.savefig("$(outpath)/ProtonDensity.png", dpi=300, transparent=true)

for idx = 1 : nFrame
    img = imgs[idx,:,:]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/recon_$(idx).png", dpi=300, transparent=true)
end



