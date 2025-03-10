using MAT
using PyPlot
using Statistics

T = Float64;

path = "$(@__DIR__)/Figures/Fig6/";
out_path     = "$(path)/out" ; if ispath(out_path) == false mkpath(out_path) end


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

zoom = [230:330, 70:170]

for idx in [1,2,3,4,5,6]
# idx = 1
subfolder = subfolders[idx];
out_path  = "$(path)/data/out/$(subfolder)" ; if ispath(out_path) == false mkpath(out_path) end

data = datas[idx]
data_file = "$(path)/data/$(data)"
@info "data file: $(data_file)"
data = matread(data_file);
imgs          = data["imgs"];
labels        = data["labels"];

figure_width       = 5/2.53999863
figure_height      = 5/2.53999863
vmaxp              = 99.9;
vminp              = 0;
color_facecolor    = "#ffffff";
dpi                = 900;

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
    ax.imshow(img[zoom[1], zoom[2]], cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(out_path)/zoom_$(label).png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.)
end
end

