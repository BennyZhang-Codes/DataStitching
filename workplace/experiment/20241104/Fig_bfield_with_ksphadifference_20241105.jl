using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

idx = 6;

path         = "E:/pulseq/20241104_ABDL/"
seqs = ["7T_1mm-200-r4_max51-fa90.seq",
        "7T_1mm-200-r3_max51-fa90.seq",
        "7T_1mm-200-r2_max51-fa90.seq",
        "7T_0.71mm-280-r2_max51-fa90.seq",
        "7T_0.6mm-332-r3_max51-fa90.seq",
        "7T_0.5mm-400-r4_max51-fa90.seq"]

dfcs = ["1p0_200_r4", "1p0_200_r3", "1p0_200_r2", "0p71_280_r2", "0p6_332_r3", "0p5_400_r4"]

seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin(seqs[idx], f)][1]
dfc_file     = "$(path)/dfc/" * [f for f in readdir("$(path)/dfc") if occursin(dfcs[idx], f)][1]


outpath = "$(path)/out/Fig"; if ispath(outpath) == false mkpath(outpath) end

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
figure_width       = 7/2.54
figure_height      = 5/2.54
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
color_kspha_diff   = "#555555"
color_bfield_diff  = "#cccccc"



fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
terms = [5,6,7,8,9]

twinxs = [axs[i].twinx() for i in 1:5]

for ax in axs[1:4, :]
    ax.tick_params(axis="both", labelbottom=false)
end
for idx in eachindex(axs)
    ax = axs[idx]
    term = terms[idx] - 1
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
bfield_vmaxs = [0.05, 55, 55, 0.2, 2, 2, 0.4, 1, 1];
kspha_vmaxs  = 1.5 .* vec(maximum(abs.(ksphaStitched .- ksphaStandard), dims=1));



for row = 1:5
    for col = 1:1
        term = terms[row]
        ax = axs[row, col]
        ax_twin = twinxs[row, col]
        kspha_diff = ksphaStitched[:, term] - ksphaStandard[:, term];
        bfield_Stitched = bfieldStitched[:, term];
        bfield_Standard = bfieldStandard[:, term];
        bfield_diff = bfieldStitched[:, term] - bfieldStandard[:, term];
        ax.set_ylim(-bfield_vmaxs[term], bfield_vmaxs[term])
        ax_twin.set_ylim(-kspha_vmaxs[term], kspha_vmaxs[term])
        ax.set_xlim(t[1], t[end])

        # line2, = ax.plot(t, bfield_diff, color=color_bfield_diff, linewidth=linewidth, label="Field difference")
        if row == 1
            line2, = ax.plot(t, bfield_Standard, color=color_bfield_diff, linewidth=linewidth, label="Standard (mT/m²)")
            line1, = ax.plot(t, bfield_Stitched, color=sh.dict["h$(term-1)"].color, linewidth=linewidth, label="$(bfield_ylabels[term])")
            line3,  = ax_twin.plot(t, kspha_diff, color=color_kspha_diff, linewidth=linewidth, label="Coefficent difference (rad/m²)")
        else
            line2, = ax.plot(t, bfield_Standard, color=color_bfield_diff, linewidth=linewidth, label="")
            line1, = ax.plot(t, bfield_Stitched, color=sh.dict["h$(term-1)"].color, linewidth=linewidth, label="$(bfield_ylabels[term])")
            line3,  = ax_twin.plot(t, kspha_diff, color=color_kspha_diff, linewidth=linewidth, label="")
        end
            # ax_twin.set_ylabel("($(units[term]))",rotation=0, 
        #     ha="left", va="center", x=0, y=0.5,
        #     fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        # ax.set_ylabel("$(bfield_ylabels[term])",rotation=0,    # \n(m$(sh.dict["h$(term-1)"].unit))
            # ha="right", va="center", x=0, y=0.5,
            # fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.legend(handles=[line1,line2], fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
            loc="upper left", bbox_to_anchor=(0,1.2),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
        ax_twin.legend(handles=[line3], fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
            loc="upper right", bbox_to_anchor=(1,1.2),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
    end
end
axs[5, 1].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)
fig.align_ylabels()
# fig.subplots_adjust(left=0.1, right=0.95, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)
fig.tight_layout(pad=0)
fig.savefig("$(outpath)/Fig_kspha_diff_$(dfcs[idx]).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_kspha_diff_$(dfcs[idx]).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)


figure_width       = 5/2.54
figure_height      = 5/2.54

figure_width       = 8/2.54
figure_height      = 8/2.54
fontsize_legend    = 9
fontsize_label     = 10
fontsize_ticklabel = 7
linewidth          = 0.75

fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
ax = axs[1,1]

vmax_k = 1.1 * maximum([maximum(abs.(ksphaStitched[:,2:3])), maximum(abs.(ksphaStandard[:,2:3]))])

ax.tick_params(axis="both", length=ticklength, width=linewidth, 
color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_color(color_label)
    ax.spines[spine].set_visible(false)
end
ax.set_facecolor(color_facecolor)

ax.plot(ksphaStitched[:,2], ksphaStitched[:,3], color="C0", linewidth=linewidth, label="Stitched")
ax.plot(ksphaStandard[:,2], ksphaStandard[:,3], color="C1", linewidth=linewidth, label="Standard")


ax.set_aspect(1)

ax.set_xlabel(L"k_x \ (rad/m)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)
ax.set_ylabel(L"k_y \ (rad/m)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)

ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
    loc="upper left", bbox_to_anchor=(0,1.05), frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)

fig.tight_layout(pad=0)
fig.savefig("$(outpath)/Fig_kspha_diff_xy_$(dfcs[idx])_8x8.png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
fig.savefig("$(outpath)/Fig_kspha_diff_xy_$(dfcs[idx])_8x8.svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)

