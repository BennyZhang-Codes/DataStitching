include("$(@__DIR__)/Figures/plt_preset.jl");
using MAT
import Statistics: quantile

outpath = "$(@__DIR__)/Figures/Fig6/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

###### load images reconstructed by different considerations of the extended signal model 
solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;
snr=10;
matfile  = "R30_snr$(snr)_$(solver)_$(iter)_$(regularization)_$(λ)"

matdatas = MAT.matread("$(outpath)/$(matfile).mat")
imgs     = matdatas["imgs"];
B0map    = matdatas["B0map"];
headmask = matdatas["headmask"];
x_ref    = matdatas["x_ref"];
nFrame, nX, nY = size(imgs)


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


###### colorbar setting for the ρ map of the phantom
ρ_values = sort(unique(x_ref))
cm_gray = plt.cm.get_cmap("gray")
cm_gray_ref = cm_gray.from_list("gray_ref",  cm_gray(ρ_values), length(ρ_values))
norm = mcolors.BoundaryNorm([-1; ρ_values] , cm_gray_ref.N)


function GT_dB0(ax)
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
end

function GT_ρ(ax)
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
end


###### plotting settings
width_ratios = [0.9, 0.5, 1, 1, 1, 1]
height = 5
width = height/2*sum(width_ratios)
@info "width = $(width), height = $(height)"

matplotlib.rc("mathtext", default="regular");
matplotlib.rc("figure", dpi=200);
matplotlib.rc("font", family="Arial");
matplotlib.rcParams["mathtext.default"];
figure_width       = width/2.53999863
figure_height      = height/2.53999863
vmaxp              = 99;
vminp              = 1;
cmap               = "gray";
fontsize_legend    = 7;
fontsize_label     = 7;
fontsize_subfigure = 9;
fontsize_ticklabel = 5;

linewidth          = 0.5;
ticklength         = 1.5;

pad_label          = 2;
pad_label_ylabel   = 4;
pad_labeltick      = 2;

color_facecolor    = "#ffffff";
color_label        = "#000000";
color_subfigure    = "#ffffff";

fig, axs = plt.subplots(nrows=2, ncols=6, figsize=(figure_width, figure_height), facecolor=color_facecolor, width_ratios=width_ratios)
ylabels = [L"w/o \ ΔB_{0}", L"w/ \ ΔB_{0}"]
titles  = ["Nominal", "Stitched 110", "Stitched 111", "Standard 111"]

idx_img = [5 6 7 8; 1 2 3 4]
for row = 1 : 2
    for col = 1 : 6
        ax = axs[row, col]
        ax.set_facecolor(color_facecolor)
        ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
        for spine in ax.spines  # "left", "right", "bottom", "top"
            ax.spines[spine].set_visible(false)
        end
        if col >= 3 
            idx = idx_img[row, col-2]
            img = imgs[idx,:,:]
            vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
            # ax.text(0.02, 0.98, "$(Char(96+idx))", fontsize=fontsize_subfigure, color=color_subfigure, transform=ax.transAxes, ha="left", va="top")
            ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
            if row == 1 ax.set_title(titles[col-2], fontsize=fontsize_label, color=color_label) end
            if col == 3 ax.set_ylabel(ylabels[row], fontsize=fontsize_label, color=color_label, labelpad=pad_label_ylabel, rotation=90, ha="center", va="center", x=0.5, y=0.5) end    
        end 

    end
end
GT_ρ(axs[1,1])
GT_dB0(axs[2,1])

fig.text(0, 1, "(a)", ha="left", va="bottom", fontsize=fontsize_subfigure, color=color_label)
fig.text(0.95*sum(width_ratios[1:2]) / sum(width_ratios) + 0, 1, "(b)", ha="right", va="bottom", fontsize=fontsize_subfigure, color=color_label)

fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975, wspace=0, hspace=0)
fig.savefig("$(outpath)/Fig6_GT_$(matfile).png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig6_GT_$(matfile).svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)

