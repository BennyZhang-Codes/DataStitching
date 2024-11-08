using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

T = Float64

path         = "E:/pulseq/20241104_ABDL/"

nameStandard = "r4_0p5_standard_rep1_wB0_standard_111_cgnr_L2_0.001_20"
nameStitched = "r4_0p5_standard_rep1_wB0_stitched_111_cgnr_L2_0.001_20"

mat_file_Stitched = "$(path)/out/" * [f for f in readdir("$(path)/out") if occursin("$(nameStitched).mat", f)][1]
mat_file_Standard = "$(path)/out/" * [f for f in readdir("$(path)/out") if occursin("$(nameStandard).mat", f)][1]
outpath = "$(path)/out/r4_0p5_standard_rep1_111_cgnr_L2_0.001_20"; if ispath(outpath) == false mkpath(outpath) end

labels = [matread(mat_file_Stitched)["label"], matread(mat_file_Standard)["label"]];


nX, nY = size(matread(mat_file_Stitched)["img"]);
imgs = Array{Complex{T},3}(undef, 2, nX, nY);
imgs[1,:,:] = matread(mat_file_Stitched)["img"];
imgs[2,:,:] = matread(mat_file_Standard)["img"];



nImg, _,_ = size(imgs);


figure_width       = 5/2.54
figure_height      = 5/2.54
vmaxp              = 99;
vminp              = 1;
color_facecolor    = "#ffffff";


vmin, vmax = quantile(abs.(imgs)[:], 0), quantile(abs.(imgs)[:], 0.999)
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
    ax.imshow(img, cmap="gray", vmin=vmin, vmax=vmax)
    fig.tight_layout(pad=0)
    fig.savefig("$(outpath)/$(label).png", dpi=300, transparent=true)
end



# difference between the stitched 111 and standard 111
scale = 10;
for scale = [5, 10, 20]
diff = abs.(abs.(imgs[1,:,:]) - abs.(imgs[2,:,:]));

vmin, vmax = quantile(abs.(imgs)[:], 0), quantile(abs.(imgs)[:], 0.999)

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1, 1]
ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end
ax.imshow(diff*scale, cmap="gray", vmin=vmin, vmax=vmax)
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/wB0_diff_111_scale$(scale).png", dpi=300, transparent=true)
end





