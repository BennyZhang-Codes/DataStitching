using PyPlot, MAT, LinearAlgebra
using HighOrderMRI
lines = PyPlot.matplotlib.lines

outpath = "$(@__DIR__)/Figures/Fig8/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

sh = SphericalHarmonics();

path = "E:/pulseq/20250222_LYQ"
seq_file_1p0 = "$(path)/seq/7T_1p0_200_r4_max51_fa90.seq"
seq_file_0p5 = "$(path)/seq/7T_0p5_400_r4_max51_fa90.seq"
dfc_file_1p0 = "$(path)/dfc_2nd/7T_2nd_tt_H15_1p0_200_r4.mat"
dfc_file_0p5 = "$(path)/dfc_2nd/7T_2nd_tt_H15_0p5_400_r4.mat"

# seq
seq_1p0 = read_seq(seq_file_1p0);
seq_0p5 = read_seq(seq_file_0p5);

# dfc
nSample_1p0, nTerm_1p0 = size(matread(dfc_file_1p0)["bfieldStitched"]);  
nSample_1p0 = Int64(round(nSample_1p0, digits=-2));
dt_1p0             = matread(dfc_file_1p0)["dt"];  # [s]
ksphaStitched_1p0  = matread(dfc_file_1p0)["ksphaStitched"][45:nSample_1p0+44, 1:nTerm_1p0]; # rad, rad/m, rad/m²
ksphaStandard_1p0  = matread(dfc_file_1p0)["ksphaStandard"][45:nSample_1p0+44, 1:nTerm_1p0]; 
bfieldStitched_1p0 = matread(dfc_file_1p0)["bfieldStitched"][1:nSample_1p0, 1:nTerm_1p0]; 
bfieldStandard_1p0 = matread(dfc_file_1p0)["bfieldStandard"][1:nSample_1p0, 1:nTerm_1p0]; # T, T/m, T/m²
ntStitched_1p0     = matread(dfc_file_1p0)["nSampleAllSegStitched"]; 
t_1p0 = collect(0:nSample_1p0-1) .* dt_1p0 .* 1e3;  # convert to ms  

a = Int64.(vec(ntStitched_1p0))
e = cumsum(a)
s = e .- a
s[1] = 1
Segment_range_1p0 = diag([st:en for st in s, en in e])


nSample_0p5, nTerm_0p5 = size(matread(dfc_file_0p5)["bfieldStitched"]);   
nSample_0p5 = Int64(round(nSample_0p5, digits=-2))
dt_0p5             = matread(dfc_file_0p5)["dt"];  # [s]
ksphaStitched_0p5  = matread(dfc_file_0p5)["ksphaStitched"][45:nSample_0p5+44, 1:nTerm_0p5]; # rad, rad/m, rad/m²
ksphaStandard_0p5  = matread(dfc_file_0p5)["ksphaStandard"][45:nSample_0p5+44, 1:nTerm_0p5]; 
bfieldStitched_0p5 = matread(dfc_file_0p5)["bfieldStitched"][1:nSample_0p5, 1:nTerm_0p5]; 
bfieldStandard_0p5 = matread(dfc_file_0p5)["bfieldStandard"][1:nSample_0p5, 1:nTerm_0p5]; # T, T/m, T/m²
ntStitched_0p5     = matread(dfc_file_0p5)["nSampleAllSegStitched"]; 
t_0p5 = collect(0:nSample_0p5-1) .* dt_0p5 .* 1e3;  # convert to ms  

a = Int64.(vec(ntStitched_0p5))
e = cumsum(a)
s = e .- a
s[1] = 1
Segment_range_0p5 = diag([st:en for st in s, en in e])

# nSegment_1p0 = length(Segment_range_1p0)
# nSegment_0p5 = length(Segment_range_0p5)
# t_1p0 = collect(range(1, nSample_1p0)) .* dt_1p0 .* 1e3;  # convert to ms  
# t_0p5 = collect(range(1, nSample_0p5)) .* dt_0p5 .* 1e3;  # convert to ms  



matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 17/2.53999863
figure_height      = 8/2.53999863
linewidth          = 0.5
linewidth_plot     = 0.3
linewidth_legend   = 0.8
ticklength         = 1.5

fontsize_subfigure = 9
fontsize_label     = 7
fontsize_legend    = 7
fontsize_ticklabel = 6

pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#ffffff"
color_label        = "#000000"
color_difference   = "#333333"
color_segment      = ["C0", "C1", "C2", "C3"]

fig, axs = plt.subplots(nrows=7, ncols=2, figsize=(figure_width, figure_height), facecolor=color_facecolor)
terms = [1, 4, 5, 6, 7, 8, 9]

for ax in axs[1:6, :]
    ax.tick_params(axis="both", labelbottom=false)
end
axs[7,1].set_xticks([0, 5, 10, 15, 20, 25, 30])
axs[7,2].set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80])
for row = 1:7
    for col = 1:2
        ax = axs[row, col]
        term = terms[row]
        ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
            color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
        # ax.tick_params(axis="y", color=sh.dict["h$(term)"].color)
        for spine in ax.spines  # "left", "right", "bottom", "top"
            ax.spines[spine].set_linewidth(linewidth)
        end
        # ax.spines["left"].set_color(sh.dict["h$(term)"].color)
        for spine in ["right", "bottom", "top"]
            ax.spines[spine].set_visible(false)
        end
        ax.set_facecolor(color_facecolor)
    end
end
ylabels = [L"G_{0}", L"G_{x}", L"G_{y}", L"G_{z}", L"G_{xy}", L"G_{zy}", L"G_{2z^2-(x^2+y^2)}", L"G_{xz}", L"G_{x^2-y^2}"]
ylabels = [L"b_{0}", L"b_{x}", L"b_{y}", L"b_{z}", L"b_{xy}", L"b_{zy}", L"b_{2z^2-(x^2+y^2)}", L"b_{xz}", L"b_{x^2-y^2}"]
units   = [L"mT", L"mT/m", L"mT/m", L"mT/m", L"mT/m^2", L"mT/m^2", L"mT/m^2", L"mT/m^2", L"mT/m^2"]


bfield_vmaxs = [0.05, 55, 55, 0.2, 2, 2, 0.5, 1.5, 1];
t_adc = [t_1p0, t_0p5];
samples_Stitched = [bfieldStitched_1p0, bfieldStitched_0p5];
samples_Standard = [bfieldStandard_1p0, bfieldStandard_0p5];
ncols = ones(14)*2;
for row = 1:7
    for col = 1:2
        ax = axs[row, col]
        term = terms[row]
        vmax = bfield_vmaxs[term]
        t = t_adc[col]
        sampleStitched = samples_Stitched[col]
        sampleStandard = samples_Standard[col]
        ncol = ncols[row]
        ax.set_ylim(-vmax, vmax)
        ax.set_xlim(t[1], t[end])

        ax.plot(t, sampleStitched[:, term] - sampleStandard[:, term], color=color_difference, linewidth=linewidth_plot)
        ax.plot(t, sampleStitched[:, term], color=sh.dict["h$(term-1)"].color, linewidth=linewidth_plot)
        ax.yaxis.set_major_locator(plt.MultipleLocator(vmax))
        ax.set_ylabel("$(ylabels[term])\n[$(units[term])]",rotation=0, 
            ha="right", va="center", x=0, y=0.5,
            fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.legend(
            [lines.Line2D([], [], color=sh.dict["h$(term)"].color, linewidth=linewidth_legend), 
                lines.Line2D([], [], color=color_difference, linewidth=linewidth_legend)], 
            ["Stitched", "Difference"], 
            handler_map=Dict(Tuple(ls) => legend_handler.HandlerTuple(ndivide=length(color_segment), pad=0)),
            fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
            loc="center left", bbox_to_anchor=(0,0.90),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
    end
end


axs[7, 1].set_xlabel("Time [ms]", fontsize=fontsize_label, color=color_label, labelpad=pad_label)
axs[7, 2].set_xlabel("Time [ms]", fontsize=fontsize_label, color=color_label, labelpad=pad_label)

fig.text(0.0, 1, "(a)", ha="left", va="bottom", fontsize=fontsize_subfigure, color=color_label)
fig.text(0.5, 1, "(b)", ha="left", va="bottom", fontsize=fontsize_subfigure, color=color_label)

fig.align_ylabels()
fig.tight_layout(pad=0, h_pad=-0.2, w_pad=0.5)

fig.savefig("$(outpath)/Fig8.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig8.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(@__DIR__)/Figures/tiff/Fig8.tiff", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)

