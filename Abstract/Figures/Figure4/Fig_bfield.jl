using PyPlot
using KomaHighOrder

sh = SphericalHarmonics()

using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

path         = "E:/pulseq/20241010_skope_fa90/invivo"
r=4
seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("r$(r)", f)][1]
dfc_file     = "$(path)/dfc/" * [f for f in readdir("$(path)/dfc") if occursin("r$(r).mat", f)][1]

outpath = "workplace/Abstract/Figures/Figure4/out"; if ispath(outpath) == false mkpath(outpath) end

dt             = matread(dfc_file)["dt"];  # [s]
Stitched  = matread(dfc_file)["bfieldStitched"]; # rad, rad/m, rad/m²
Standard  = matread(dfc_file)["bfieldStandard"]; 
nSample, nTerm = size(matread(dfc_file)["bfieldStitched"]);   

t = collect(0:nSample-1) .* dt .* 1e3;  # convert to ms  


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
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
color_difference   = "#999999"

fig, axs = plt.subplots(nrows=9, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)

for ax in axs[1:8, :]
    ax.tick_params(axis="both", labelbottom=false)
end
for ax in axs
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_linewidth(linewidth)
        # if spine in ["right", "top"]
        #     ax.spines[spine].set_visible(false)
        # end
    end
    ax.set_facecolor(color_facecolor)
end
ylabels = [L"G_{0}", L"G_{x}", L"G_{y}", L"G_{z}", L"G_{xy}", L"G_{zy}", L"G_{2z^2-(x^2+y^2)}", L"G_{xz}", L"G_{x^2-y^2}"]
# ylabels = [L"k_{0}", L"k_{x}", L"k_{y}", L"k_{z}", L"k_{xy}", L"k_{zy}", L"k_{2z^2-(x^2+y^2)}", L"k_{xz}", L"k_{x^2-y^2}"]

vmaxs = [0.05, 50, 50, 0.2, 2, 2, 0.4, 1, 1];

for row = 1:9
    for col = 1:1
        ax = axs[row, col]
        term = row-1
        vmax = vmaxs[row]
        sampleStitched = Stitched[:, row]
        sampleStandard = Standard[:, row]
        ax.set_ylim(-vmax, vmax)
        ax.set_xlim(t[1], t[end])

        line1, = ax.plot(t, (sampleStitched .- sampleStandard), color=color_difference, linewidth=linewidth, label="Difference")
        line2, = ax.plot(t, sampleStitched, color=sh.dict["h$(term)"].color, linewidth=linewidth, label="Stitched measurement")
        # ax.yaxis.set_major_locator(plt.MultipleLocator(vmax))
        ax.set_ylabel("$(ylabels[term+1])\n(m$(sh.dict["h$(term)"].unit))",rotation=0, 
            ha="right", va="center", x=0, y=0.5,
            fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.legend(handles=[line2, line1], fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
            loc="upper left", bbox_to_anchor=(0,1.15),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
    end
end
axs[9, 1].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)
fig.align_ylabels()
# fig.subplots_adjust(left=0.1, right=0.95, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)
fig.tight_layout(pad=0, h_pad=0.1, w_pad=0.3)
fig.savefig("$(outpath)/Fig_bfield_r$(r).png", dpi=300, bbox_inches="tight", transparent=true)
fig.savefig("$(outpath)/Fig_bfield_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true)


