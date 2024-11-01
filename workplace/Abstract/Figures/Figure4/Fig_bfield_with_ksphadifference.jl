using PyPlot
using KomaHighOrder

sh = SphericalHarmonics()

using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

path         = "E:/pulseq/20241010_skope_fa90/invivo"
r=2
seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("r$(r)", f)][1]
dfc_file     = "$(path)/dfc/" * [f for f in readdir("$(path)/dfc") if occursin("r$(r).mat", f)][1]

outpath = "workplace/Abstract/Figures/Figure4/out"; if ispath(outpath) == false mkpath(outpath) end

dt             = matread(dfc_file)["dt"];  # [s]
ksphaStitched  = matread(dfc_file)["ksphaStitched"]; # rad, rad/m, rad/m²
ksphaStandard  = matread(dfc_file)["ksphaStandard"]; 
bfieldStitched = matread(dfc_file)["bfieldStitched"]; 
bfieldStandard = matread(dfc_file)["bfieldStandard"]; # T, T/m, T/m²

nSample, nTerm = size(matread(dfc_file)["bfieldStitched"]);   
ksphaStandard = ksphaStandard[end-nSample+1:end, 1:nTerm];
ksphaStitched = ksphaStitched[end-nSample+1:end, 1:nTerm];


t = collect(0:nSample-1) .* dt .* 1e3;  # convert to ms  


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 9/2.54
figure_height      = 9/2.54
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
color_kspha_diff   = "#777777"
color_bfield_diff  = "#bbbbbb"

fig, axs = plt.subplots(nrows=9, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
twinxs = [axs[i].twinx() for i in 1:9]

for ax in axs[1:8, :]
    ax.tick_params(axis="both", labelbottom=false)
end
for idx in eachindex(axs)
    ax = axs[idx]
    term = idx-1
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    ax.tick_params(axis="y", color=sh.dict["h$(term)"].color)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_linewidth(linewidth)
    end
    ax.spines["left"].set_color(sh.dict["h$(term)"].color)
    for spine in ["right", "bottom", "top"]
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
end

for ax in twinxs
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_linewidth(linewidth)
    end
    ax.tick_params(axis="y", color=color_kspha_diff)
    ax.spines["right"].set_color(color_kspha_diff)
    for spine in ["left", "bottom", "top"]
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
end
bfield_ylabels = [L"b_{0}", L"b_{x}", L"b_{y}", L"b_{z}", L"b_{xy}", L"b_{zy}", L"b_{2z^2-(x^2+y^2)}", L"b_{xz}", L"b_{x^2-y^2}"]
kspha_ylabels  = [L"k_{0}", L"k_{x}", L"k_{y}", L"k_{z}", L"k_{xy}", L"k_{zy}", L"k_{2z^2-(x^2+y^2)}", L"k_{xz}", L"k_{x^2-y^2}"]
units = ["rad", "rad/m", "rad/m", "rad/m", "rad/m²", "rad/m²", "rad/m²", "rad/m²", "rad/m²"]
bfield_vmaxs = [0.05, 50, 50, 0.2, 2, 2, 0.4, 1, 1];
# bfield_vmaxs = 1.2 .* vec(maximum(abs.(bfieldStitched), dims=1));
kspha_vmaxs  = 1.5 .* vec(maximum(abs.(ksphaStitched .- ksphaStandard), dims=1));

for row = 1:9
    for col = 1:1
        ax = axs[row, col]
        ax_twin = twinxs[row, col]
        term = row-1
        kspha_diff = ksphaStitched[:, row] - ksphaStandard[:, row];
        bfield_Stitched = bfieldStitched[:, row];
        bfield_diff = bfieldStitched[:, row] - bfieldStandard[:, row];
        ax.set_ylim(-bfield_vmaxs[row], bfield_vmaxs[row])
        ax_twin.set_ylim(-kspha_vmaxs[row], kspha_vmaxs[row])
        ax.set_xlim(t[1], t[end])

        # line2, = ax.plot(t, bfield_diff, color=color_bfield_diff, linewidth=linewidth, label="Field difference")
        line1, = ax.plot(t, bfield_Stitched, color=sh.dict["h$(term)"].color, linewidth=linewidth, label="Stitching")
        line3,  = ax_twin.plot(t, kspha_diff, color=color_kspha_diff, linewidth=linewidth, label="Coefficent difference")
        ax_twin.set_ylabel("$(kspha_ylabels[term+1])\n($(units[row]))",rotation=0, 
            ha="left", va="center", x=0, y=0.5,
            fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.set_ylabel("$(bfield_ylabels[term+1])\n(m$(sh.dict["h$(term)"].unit))",rotation=0, 
            ha="right", va="center", x=0, y=0.5,
            fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.legend(handles=[line1,], fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
            loc="upper left", bbox_to_anchor=(0,1.3),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
        ax_twin.legend(handles=[line3], fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
            loc="upper right", bbox_to_anchor=(1,1.3),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
    end
end
axs[9, 1].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)
fig.align_ylabels()
# fig.subplots_adjust(left=0.1, right=0.95, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/Fig_kspha_diff_r$(r).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_kspha_diff_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)


