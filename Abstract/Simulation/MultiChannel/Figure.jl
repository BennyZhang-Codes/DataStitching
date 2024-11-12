include("Figures/Fig_preset.jl");
using MAT
import Statistics: quantile


matfile = "fully_snr10_admm_20_TV_0.0001"
path = "$(@__DIR__)/Abstract/Simulation/MultiChannel/out"; if ispath(path) == false mkpath(path) end     # output directory
outpath = "$(outpath)/$(matfile)"; if ispath(outpath) == false mkpath(outpath) end     # output directory)";
###### load images reconstructed by different considerations of the extended signal model 
imgs = MAT.matread("$(path)/$(matfile).mat")["imgs"];
labels = MAT.matread("$(path)/$(matfile).mat")["titles"];

# csm      = results["csm"];
# B0map    = results["B0map"];
# headmask = results["headmask"];
# x_ref    = results["x_ref"];
# imgs     = results["imgs"];
# signal   = results["signal"];
# labels   = results["labels"];

matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=100)
matplotlib.rc("font", family="Arial")    # "Times New Roman"
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

vmaxp              = 99.0;
vminp              = 1;

#############################################################################
# plot signal in each channel
#############################################################################
nImg, nX, nY = size(imgs)
vmin, vmax = quantile(imgs[4:4,:,:][:], vminp/100), quantile(imgs[4:4,:,:][:], vmaxp/100)
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
    fig.tight_layout(pad=0)
    # fig.savefig("$(outpath)/img_$(idx).png", dpi=300, transparent=true)
end


for scale = [1,2,3]
    diff = abs.(abs.(imgs[7,:,:]) - abs.(imgs[8,:,:]));
    
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    vmin, vmax = quantile(imgs[4:4,:,:][:], vminp/100), quantile(imgs[4:4,:,:][:], vmaxp/100)
    ax.imshow(diff*scale, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/wB0_diff_111_scale$(scale).png", dpi=300, transparent=true)
end

#############################################################################
# plot signal in each channel
#############################################################################
figure_width       = 3.5/2.54
figure_height      = 1.8/2.54

nSample, nCha = size(data);
for cha = 1 : nCha
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor)

    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
    ax.plot(abs.(data[:, cha]), linewidth=0.5, color="C$(cha%9)")
    fig.tight_layout(pad=0)

    # fig.savefig("$(outpath)/Fig_signal_r$(r)_cha$(cha).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
    # fig.savefig("$(outpath)/Fig_signal_r$(r)_cha$(cha).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
end
