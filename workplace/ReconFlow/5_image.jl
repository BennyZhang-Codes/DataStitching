using MAT
using PyPlot
using Statistics

T = Float64;

path         = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/"
out_path     = "$(path)/out" ; if ispath(out_path) == false mkpath(out_path) end


datas = ["1p0_200_r4.mat",
         "1p0_200_r3.mat",
         "1p0_200_r2.mat",
         "0p71_280_r2.mat",
         "0p6_332_r3.mat",
         "0p5_400_r4.mat"]

for idx in [1,2,3,4,5,6]
# idx = 1

subfolder = datas[idx][1:end-4];
out_path  = "$(path)/out/$(subfolder)" ; if ispath(out_path) == false mkpath(out_path) end

data = datas[idx]
data_file = "$(path)/data/$(data)"
@info "data file: $(data_file)"
data = matread(data_file);

dt            = data["seq_dt"];

tauStitched   = data["dfc_tauStitched"];
tauStandard   = data["dfc_tauStandard"];
tauNominal    = data["tauNominal"];

recon         = data["recon"];
imgStitched   = recon["imgStitched"];
imgStandard   = recon["imgStandard"];
imgNominal    = recon["imgNominal"];

labelStitched = recon["labelStitched"];
labelStandard = recon["labelStandard"];
labelNominal  = recon["labelNominal"];


figure_width       = 5/2.53999863
figure_height      = 5/2.53999863
vmaxp              = 99.9;
vminp              = 0;
color_facecolor    = "#ffffff";
dpi                = 900;


imgs   = imgStitched
labels = labelStitched
nImg   = length(labels)
vmin, vmax = quantile(abs.(imgs)[:], vminp/100), quantile(abs.(imgs)[:], vmaxp/100)
for idx = 1 : nImg
    img = abs.(imgs)[:,:, idx]
    label = labels[idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(out_path)/$(label).png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
end
fig = plt_images(  abs.(imgs); dim=3, height=5, width=5, vmaxp=vmaxp, vminp=vminp)
fig.savefig("$(out_path)/Stitched_mag.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
fig = plt_images(angle.(imgs); dim=3, height=5, width=5, vmax=π, vmin=-π)
fig.savefig("$(out_path)/Stitched_pha.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)


imgs   = imgStandard
labels = labelStandard
nImg   = length(labels)
vmin, vmax = quantile(abs.(imgs)[:], vminp/100), quantile(abs.(imgs)[:], vmaxp/100)
for idx = 1 : nImg
    img = abs.(imgs)[:,:, idx]
    label = labels[idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(out_path)/$(label).png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
end
fig = plt_images(  abs.(imgs); dim=3, height=5, width=5, vmaxp=vmaxp, vminp=vminp)
fig.savefig("$(out_path)/Standard_mag.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
fig = plt_images(angle.(imgs); dim=3, height=5, width=5, vmax=π, vmin=-π)
fig.savefig("$(out_path)/Standard_pha.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)


imgs   = imgNominal
labels = labelNominal
nImg   = length(labels)
vmin, vmax = quantile(abs.(imgs)[:], vminp/100), quantile(abs.(imgs)[:], vmaxp/100)
for idx = 1 : nImg
    img = abs.(imgs)[:,:, idx]
    label = labels[idx]
    fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
    ax = axs[1, 1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(out_path)/$(label).png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
end

end