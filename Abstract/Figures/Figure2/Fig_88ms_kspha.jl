using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

r=1
# gre_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(r"^syn.*gre.*mat$", f)][1]
seq_file = "$(@__DIR__)/demo/SingleChannel/1mm_R1.seq"        # *.seq file is the pulseq's sequence file
dfc_file = "$(@__DIR__)/demo/SingleChannel/1mm_R1.mat"   # *.mat file contains the dynamic field data from both stitching method and the standard method.
outpath = "workplace/Abstract/Figures/Figure2/out"; if ispath(outpath) == false mkpath(outpath) end

print("$(@__DIR__)")
########################################################################
# 1. load the *.seq file and extract some infomations
########################################################################
seq = read_seq(seq_file)[4:end]; 
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
t_gradFreeDelay = seq.DEF["skope_gradFreeDelay"]; # [s]

########################################################################
# 2. load the *.mat file and extract the DFC data
########################################################################
dt            = matread(dfc_file)["dt"];  # [s]
bfieldStitched = matread(dfc_file)["skopeStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
bfieldStandard = matread(dfc_file)["skopeStandard"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
nSampleAllSeg = matread(dfc_file)["ntStitched"]; 
nGradPoint, nTerm = size(matread(dfc_file)["skopeStitched"]);   

t = dt * (nGradPoint-1);
GR_bfieldStitched = reshape([KomaMRIBase.Grad(bfieldStitched[idx,:], t, dt/2, dt/2, 0) for idx=1:9], :, 1);
GR_bfieldStandard = reshape([KomaMRIBase.Grad(bfieldStandard[idx,:], t, dt/2, dt/2, 0) for idx=1:9], :, 1);

hoseqStitched = HO_Sequence(seq);             # hoseq, defined in KomaHighOrder.jl
hoseqStandard = HO_Sequence(seq);             # hoseq, defined in KomaHighOrder.jl
hoseqStitched.GR_dfc[2:4, :] = hoseqStitched.SEQ.GR;  # copy the 1st-order gradient data from the seq object to the hoseq object
hoseqStandard.GR_dfc[2:4, :] = hoseqStandard.SEQ.GR;  # copy the 1st-order gradient data from the seq object to the hoseq object
hoseqStitched.GR_dfc[:,5] = GR_bfieldStitched;           # "5" is the index of the readout block in the spiral sequence
hoseqStandard.GR_dfc[:,5] = GR_bfieldStandard;           # "5" is the index of the readout block in the spiral sequence

# plot_seq(hoseqStitched)
# plot_seq(hoseqStandard)
# finally, hoseq* contains both the nominal trajectory and the measured trajectory (up to 2nd-order)

_, _, _, ksphaStitched = get_kspace(hoseqStitched; Δt=1)
_, _, _, ksphaStandard = get_kspace(hoseqStandard; Δt=1)

nSample, nTerm = size(ksphaStitched);   

Stitched = ksphaStandard * 2π;
t = collect(0:nSample-1) .* dt .* 1e3;  # convert to ms  


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
# matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 4.5/2.54
figure_height      = 3/2.54
linewidth          = 0.5
ticklength         = 1.5
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
fontsize_subfigure = 8
pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#ffffff"
color_label        = "#000000"
color_difference   = "#999999"

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false, sharex=true)
bfield_ylabels = [L"k_{0}", L"k_{x}", L"k_{y}", L"k_{z}", L"k_{xy}", L"k_{zy}", L"k_{2z^2-(x^2+y^2)}", L"k_{xz}", L"k_{x^2-y^2}"]
y_labels = ["Zeroth order\n(rad)", "First order\n(rad/m)", "Second order\n(rad/m²)"]

for row = 1:3
    col = 1
    ax = axs[row, col]
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    # ax.set_ylabel(y_labels[row], rotation=0, ha="right", va="center", x=0, y=0.5, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    # ax.spines["bottom"].set_visible(true)
    ax.set_facecolor(color_facecolor)
end

order = [1,2,2,2,3,3,3,3,3];


for term = 0:0
    idx = term+1
    ax = axs[order[idx], 1]
    ymax = 1.1*maximum(abs.(Stitched[:, 1]))
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=1, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for term = 1:3
    idx = term+1
    ax = axs[order[idx], 1]
    ymax = 1.1*maximum(abs.(Stitched[:, 2:4]))
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=3, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for term = 4:8
    idx = term+1
    ax = axs[order[idx], 1]
    ymax = 1.1*maximum(abs.(Stitched[:, 5:9]))
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=5, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end


fig.align_ylabels()
fig.tight_layout(pad=0, h_pad=0)
fig.savefig("$(outpath)/Fig_kspha_standard_r$(r).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_kspha_standard_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)


Stitched = ksphaStitched * 2π;
t = collect(0:nSample-1) .* dt .* 1e3;  # convert to ms  

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false, sharex=true)
bfield_ylabels = [L"k_{0}", L"k_{x}", L"k_{y}", L"k_{z}", L"k_{xy}", L"k_{zy}", L"k_{2z^2-(x^2+y^2)}", L"k_{xz}", L"k_{x^2-y^2}"]
y_labels = ["Zeroth order\n(rad)", "First order\n(rad/m)", "Second order\n(rad/m²)"]

for row = 1:3
    col = 1
    ax = axs[row, col]
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    # ax.set_ylabel(y_labels[row], rotation=0, ha="right", va="center", x=0, y=0.5, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    # ax.spines["bottom"].set_visible(true)
    ax.set_facecolor(color_facecolor)
end

order = [1,2,2,2,3,3,3,3,3];


for term = 0:0
    idx = term+1
    ax = axs[order[idx], 1]
    ymax = 1.1*maximum(abs.(Stitched[:, 1]))
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=1, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for term = 1:3
    idx = term+1
    ax = axs[order[idx], 1]
    ymax = 1.1*maximum(abs.(Stitched[:, 2:4]))
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=3, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for term = 4:8
    idx = term+1
    ax = axs[order[idx], 1]
    ymax = 1.1*maximum(abs.(Stitched[:, 5:9]))
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=5, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end


fig.align_ylabels()
fig.tight_layout(pad=0, h_pad=0)
fig.savefig("$(outpath)/Fig_kspha_stitched_r$(r).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_kspha_stitched_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
