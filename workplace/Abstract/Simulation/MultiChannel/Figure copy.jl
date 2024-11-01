include("Figures/Fig_preset.jl");
using MAT
import Statistics: quantile

path = "workplace/Abstract/Simulation/MultiChannel/out"; if ispath(dir) == false mkpath(dir) end     # output directory


###### load images reconstructed by different considerations of the extended signal model 

outpath = "$(path)/diff"; if ispath(outpath) == false mkpath(outpath) end
imgs       = MAT.matread("$(path)/fully_wonoise_admm_20_TV_0.0001.mat")["imgs"];
imgs_noise = MAT.matread("$(path)/fully_snr10_admm_20_TV_0.0001.mat")["imgs"];
labels = MAT.matread("$(path)/fully_snr10_admm_20_TV_0.0001.mat")["titles"];


figure_width       = 5/2.54
figure_height      = 5/2.54
vmaxp              = 99;
vminp              = 1;
color_facecolor    = "#ffffff";




nImg, nX, nY = size(imgs)

vmin, vmax = quantile(abs.(imgs)[:], 0), quantile(abs.(imgs)[:], 0.999)
scale = 20;


for idx = 1 : nImg
    diff = abs.(abs.(imgs[idx,:,:]) - abs.(imgs_noise[idx,:,:]));
    label = labels[idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    vmin, vmax = quantile(img[:], vminp/100), quantile(img[:], vmaxp/100)
    ax.imshow(diff*scale, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/diff_$(idx).png", dpi=300, transparent=true)
end
