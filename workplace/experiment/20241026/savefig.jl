using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

T = Float64

path         = "E:/pulseq/20241010_skope_fa90/invivo"

name = "r4_cgnr_L2_0.01_10"

mat_file = "$(path)/out/" * [f for f in readdir("$(path)/out") if occursin("$(name).mat", f)][1]
outpath = "$(path)/out/" * name; if ispath(outpath) == false mkpath(outpath) end

labels = matread(mat_file)["labels"];
imgs   = matread(mat_file)["imgs"];
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
scale = 20;
diff = abs.(abs.(imgs[4,:,:]) - abs.(imgs[7,:,:]));

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

# difference between the stitched 111 and standard 111
scale = 20;
diff = abs.(abs.(imgs[11,:,:]) - abs.(imgs[14,:,:]));

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
fig.savefig("$(outpath)/woB0_diff_111_scale$(scale).png", dpi=300, transparent=true)




