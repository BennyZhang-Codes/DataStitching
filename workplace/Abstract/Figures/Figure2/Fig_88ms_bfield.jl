using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

r=1
# gre_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(r"^syn.*gre.*mat$", f)][1]
seq_file = "$(@__DIR__)/demo/Sequence/1mm_R1.seq"        # *.seq file is the pulseq's sequence file
dfc_file = "$(@__DIR__)/demo/DynamicFields/1mm_R1.mat"   # *.mat file contains the dynamic field data from both stitching method and the standard method.
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
dt             = matread(dfc_file)["dt"];  # [s]
bfieldStitched = matread(dfc_file)["skopeStitched"];
bfieldStandard = matread(dfc_file)["skopeStandard"];
nSampleAllSeg  = matread(dfc_file)["ntStitched"]; 
nSample, nTerm = size(matread(dfc_file)["skopeStitched"]);   

Stitched = bfieldStitched
t = collect(0:nSample-1) .* dt .* 1e3;  # convert to ms  



matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
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
bfield_ylabels = [L"b_{0}", L"b_{x}", L"b_{y}", L"b_{z}", L"b_{xy}", L"b_{zy}", L"b_{2z^2-(x^2+y^2)}", L"b_{xz}", L"b_{x^2-y^2}"]
y_labels = ["Zeroth order\n(mT)", "First order\n(mT/m)", "Second order\n(mT/m²)"]

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
    ax.set_facecolor(color_facecolor)
end

order = [1,2,2,2,3,3,3,3,3];
bfield_vmaxs = [0.07, 60, 4];


for term = 0:0
    idx = term+1
    ax = axs[order[idx], 1]
    ax.set_ylim(-bfield_vmaxs[order[idx]], bfield_vmaxs[order[idx]])
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=1, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for term = 1:3
    idx = term+1
    ax = axs[order[idx], 1]
    ax.set_ylim(-bfield_vmaxs[order[idx]], bfield_vmaxs[order[idx]])
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=3, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for term = 4:8
    idx = term+1
    ax = axs[order[idx], 1]
    ax.set_ylim(-bfield_vmaxs[order[idx]], bfield_vmaxs[order[idx]])
    ax.set_xlim(t[1], t[end])
    ax.plot(t, Stitched[:, idx], color=sh.dict["h$(term)"].color, linewidth=linewidth, label=bfield_ylabels[idx])
    # ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=5, 
    #     loc="upper left", bbox_to_anchor=(0,1.05),
    #     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end


fig.align_ylabels()
fig.tight_layout(pad=0, h_pad=0)
fig.savefig("$(outpath)/Fig_bfield_r$(r).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_bfield_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)


