using MAT
import Statistics: quantile
using PyPlot, PyCall
matplotlib.rcParams["mathtext.default"]
matplotlib.rc("mathtext", default="regular")
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
mcolors = matplotlib.colors

path = "$(@__DIR__)/Figures/Fig6/";
out_path     = "$(path)/out" ; if ispath(out_path) == false mkpath(out_path) end

# datas = [
#          "0p3_500_r30_pha_snr10_admm_20_TV_0.5.mat",
#          "0p3_500_r30_snr10_admm_20_TV_0.5.mat",
#          "0p3_500_r30_pha_snr10_cgnr_10_L2_0.0.mat",
#          "0p3_500_r30_pha_snr10_cgnr_20_L2_0.0.mat",
#          ]
# subfolders = [
#         "0p3_500_r30_pha_snr10_admm_20_TV_0p5",
#         "0p3_500_r30_snr10_admm_20_TV_0p5",
#         "0p3_500_r30_pha_snr10_cgnr_10_L2_0p0",
#         "0p3_500_r30_pha_snr10_cgnr_20_L2_0p0",
#          ]

datas = [
         "0p3_500_r30_144_ss1_snr10_cgnr_20_L2_1.0e-9.mat",
         "0p3_500_r30_144_ss3_snr10_cgnr_20_L2_1.0e-9.mat",
         "0p3_500_r30_196_ss1_snr10_cgnr_20_L2_1.0e-9.mat",
         "0p3_500_r30_196_ss3_snr10_cgnr_20_L2_1.0e-9.mat",
         "0p3_500_r30_256_ss1_snr10_cgnr_20_L2_1.0e-9.mat",
         "0p3_500_r30_256_ss3_snr10_cgnr_20_L2_1.0e-9.mat",
         ]
subfolders = [
        "0p3_500_r30_144_ss1_snr10_cgnr_20_L2_1e-9.mat",
        "0p3_500_r30_144_ss3_snr10_cgnr_20_L2_1e-9.mat",
        "0p3_500_r30_196_ss1_snr10_cgnr_20_L2_1e-9.mat",
        "0p3_500_r30_196_ss3_snr10_cgnr_20_L2_1e-9.mat",
        "0p3_500_r30_256_ss1_snr10_cgnr_20_L2_1e-9.mat",
        "0p3_500_r30_256_ss3_snr10_cgnr_20_L2_1e-9.mat",
         ]

for idx in [1, 2, 3, 4, 5, 6]
# idx = 4
subfolder = subfolders[idx];
out_path  = "$(path)/data/out/$(subfolder)" ; if ispath(out_path) == false mkpath(out_path) end

data = datas[idx]
data_file = "$(path)/data/$(data)"
@info "data file: $(data_file)"
data = matread(data_file);
imgs          = data["imgs"];
labels        = data["labels"];

x_ref         = data["x_ref"];
headmask      = data["headmask"];
B0map         = data["b0"];

# figure_width       = 5/2.53999863
# figure_height      = 5/2.53999863
fig_width          = 5
fig_height         = 5
vmaxp              = 99;
vminp              = 0;
color_facecolor    = "#ffffff";
dpi                = 900;

trajs = [data["ksphaNominal"], data["ksphaStitched"][:, 2:4], data["ksphaStandard"][:, 2:4]];
fig = plt_traj(trajs; labels=["Nominal", "Stitched", "Standard"], linewidth = 0.5,
        fontsize_label=9, fontsize_legend=7, fontsize_ticklabel=7)
fig.savefig("$(out_path)/trajs.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)

# imgs   = imgStitched
# labels = labelStitched
nImg   = length(labels)
# vmin, vmax = quantile(abs.(imgs)[:], vminp/100), quantile(abs.(imgs)[:], vmaxp/100)
vmin, vmax = 0, 1
for idx = 1 : nImg
    # img = abs.(imgs)[:,:, idx]
    label = labels[idx]
    img = HO_img_scale(x_ref, abs.(imgs[:,:,idx]))
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    fig = plt_image(abs.(img); height=fig_height, width=fig_width, vmax=vmax, vmin=vmin)
    fig.savefig("$(out_path)/mag_$(label).png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
    fig = plt_image(angle.(imgs[:,:,idx]); height=fig_height, width=fig_width, vmax=π, vmin=-π)
    fig.savefig("$(out_path)/pha_$(label).png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
end

fig = plt_images(  abs.(imgs); dim=3, height=fig_height, width=fig_width, vmaxp=vmaxp, vminp=vminp)
fig.savefig("$(out_path)/Recon_mag.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
fig = plt_images(angle.(imgs); dim=3, height=fig_height, width=fig_width, vmax=π, vmin=-π)
fig.savefig("$(out_path)/Recon_pha.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)

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
# matplotlib.rc("figure", dpi=200)
# matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 5/2.53999863
figure_height      = 5/2.53999863
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
# ax.set_title("Off-resonance", fontsize=fontsize_label, color=color_label)
divider = mpl_axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(ai, cax=cax)
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
cb.outline.set_visible(false)
cb.update_ticks()
fig.tight_layout(pad=0)
fig.savefig("$(out_path)/B0map.png", dpi=dpi, transparent=true)
# fig.savefig("$(outpath)/B0map.svg", dpi=dpi, transparent=true)

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
# ax.set_title("Proton density", fontsize=fontsize_label, color=color_label)
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
fig.savefig("$(out_path)/ProtonDensity.png", dpi=dpi, transparent=true)
# fig.savefig("$(outpath)/ProtonDensity.svg", dpi=dpi, transparent=true)
end