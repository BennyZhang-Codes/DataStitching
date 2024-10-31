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
# syn_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin("r$(r)", f)][1]
# gre_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(r"^syn.*gre.*mat$", f)][1]

dt             = matread(dfc_file)["dt"];  # [s]
Stitched  = matread(dfc_file)["ksphaStitched"]; # rad, rad/m, rad/m²
Standard  = matread(dfc_file)["ksphaStandard"]; 
nSample, nTerm = size(matread(dfc_file)["ksphaStitched"]);   

t = collect(0:nSample-1) .* dt .* 1e3;  # convert to ms  



matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 5/2.54
figure_height      = 2.8/2.54
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

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)


y_labels = [L"0^{th}\ order", L"1^{st}\ order", L"2^{nd}\ order"]

for row = 1:3
    col = 1
    ax = axs[row, col]
    # ax.set_ylabel(y_labels[row], rotation=0, ha="right", va="center", x=0, y=0.5, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
end

order = [1,2,2,2,3,3,3,3,3];
for term = 0:8
    idx = term+1
    ax = axs[order[idx], 1]
    sampleStitched = Stitched[:, idx]
    ax.plot(t, sampleStitched, color=sh.dict["h$(term)"].color, linewidth=linewidth)
end


fig.align_ylabels()
fig.tight_layout(pad=0)
fig.savefig("workplace/Abstract/Figure1/out/Fig_kspha_r$(r).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("workplace/Abstract/Figure1/out/Fig_kspha_r$(r).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)


