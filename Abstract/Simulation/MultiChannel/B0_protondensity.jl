include("Figures/Fig_preset.jl");
using MAT
import Statistics: quantile


outpath = "$(@__DIR__)/Abstract/Simulation/MultiChannel/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

###### get the reference images
simtype = SimType(B0=true, T2=false, ss=5);
location = 0.8;
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2); # decide which phantom file to use

db0_type = :quadratic;
db0_max = 100.;        
# ΔB₀ map
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (1., 1.); location=location, ss=simtype.ss, db0_type=db0_type, db0_max=db0_max);
B0map = rotl90(B0map);

x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=simtype.ss);
x_ref = rotl90(x_ref);

headmask = brain_phantom2D_reference(phantom, :headmask, (150., 150.), (1., 1.); location=location, ss=simtype.ss);
headmask = rotl90(headmask);

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




vmaxp              = 99;
vminp              = 1;

matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
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
# divider = mpl_axes_grid1.make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.02)
# # ax.set_title("Proton density", fontsize=fontsize_label, color=color_label)
# cax.set_title("[a.u.]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
# cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm_gray_ref, norm=norm), cax=cax, drawedges=true)
# cb.minorticks_off()
# cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel, length=ticklength, width=linewidth, pad=pad_labeltick)
# cb.set_ticks(ρ_values - (ρ_values - [-1; ρ_values][1:end-1]) / 2)
# cb.set_ticklabels([string(round(v, digits=2)) for v in ρ_values])
# cb.ax.tick_params(size=0)
# # cb.outline.set_visible(false)
# cb.outline.set_edgecolor(color_label)
# cb.outline.set_linewidth(linewidth)
# cb.dividers.set_color(color_label)
# cb.dividers.set_linewidth(linewidth)
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/ProtonDensity.png", dpi=900, transparent=true)
fig.savefig("$(outpath)/ProtonDensity.svg", dpi=900, transparent=true)

