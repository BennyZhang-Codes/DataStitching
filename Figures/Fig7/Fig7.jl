using PyPlot
using HighOrderMRI
using LinearAlgebra
using MAT
legend_handler = PyPlot.matplotlib.legend_handler
collections = PyPlot.matplotlib.collections
lines = PyPlot.matplotlib.lines
Patch = matplotlib.patches.Patch

sh = SphericalHarmonics()

outpath = "$(@__DIR__)/Figures/Fig7/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

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
# nSample_1p0 = Int64(round(nSample_1p0, digits=-2));
dt_1p0             = matread(dfc_file_1p0)["dt"];  # [s]
ksphaStitched_1p0  = matread(dfc_file_1p0)["ksphaStitched"][21:nSample_1p0+20, 1:nTerm_1p0]; # rad, rad/m, rad/m²
ksphaStandard_1p0  = matread(dfc_file_1p0)["ksphaStandard"][21:nSample_1p0+20, 1:nTerm_1p0]; 
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
# nSample_0p5 = Int64(round(nSample_0p5, digits=-2))
dt_0p5             = matread(dfc_file_0p5)["dt"];  # [s]
ksphaStitched_0p5  = matread(dfc_file_0p5)["ksphaStitched"][21:nSample_0p5+20, 1:nTerm_0p5]; # rad, rad/m, rad/m²
ksphaStandard_0p5  = matread(dfc_file_0p5)["ksphaStandard"][21:nSample_0p5+20, 1:nTerm_0p5]; 
bfieldStitched_0p5 = matread(dfc_file_0p5)["bfieldStitched"][1:nSample_0p5, 1:nTerm_0p5]; 
bfieldStandard_0p5 = matread(dfc_file_0p5)["bfieldStandard"][1:nSample_0p5, 1:nTerm_0p5]; # T, T/m, T/m²
ntStitched_0p5     = matread(dfc_file_0p5)["nSampleAllSegStitched"]; 
t_0p5 = collect(0:nSample_0p5-1) .* dt_0p5 .* 1e3;  # convert to ms  

a = Int64.(vec(ntStitched_0p5))
e = cumsum(a)
s = e .- a
s[1] = 1
Segment_range_0p5 = diag([st:en for st in s, en in e])

nSegment_1p0 = length(Segment_range_1p0)
nSegment_0p5 = length(Segment_range_0p5)
t_1p0 = collect(range(1, nSample_1p0)) .* dt_1p0 .* 1e3;  # convert to ms  
t_0p5 = collect(range(1, nSample_0p5)) .* dt_0p5 .* 1e3;  # convert to ms  




matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 17.5/2.53999863
figure_height      = 10/2.53999863
linewidth          = 0.5
linewidth_plot     = 0.6
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
color_segment      = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",  "#e377c2",  "#bcbd22","#17becf"]#,"#8c564b"  ,"#7f7f7f"
color_gx           = "C0"
color_gy           = "C1"

# color_segment      = ["#6699CC","#FF6666","#669966","#CC66CC","#FF9900"]#,"#8c564b"  ,"#7f7f7f"
# color_gx           = "#6699CC"
# color_gy           = "#FF6666"

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(figure_width, figure_height), facecolor=color_facecolor)
ax1, ax2, ax3, ax4, ax5, ax6 = axs 
vmax_g_1p0 = 1.1 * maximum(abs.(bfieldStitched_1p0[:, 2]))
vmax_k_1p0 = 1.1 * maximum([maximum(abs.(ksphaStitched_1p0[:,2:3])), maximum(abs.(ksphaStitched_1p0[:,2:3]))])
vmax_g_0p5 = 1.1 * maximum(abs.(bfieldStitched_0p5[:, 2]))
vmax_k_0p5 = 1.1 * maximum([maximum(abs.(ksphaStitched_0p5[:,2:3])), maximum(abs.(ksphaStitched_0p5[:,2:3]))])


ax1.set_xticks([0,5,10,15,20,25,30])
ax3.set_xticks([0,5,10,15,20,25,30])
ax2.set_xticks([0,20,40,60,80])
ax4.set_xticks([0,20,40,60,80])

ax5.set_xticks([-500, -250, 0, 250, 500])
ax5.set_yticks([-500, -250, 0, 250, 500])
ax6.set_xticks([-1000, -500, 0, 500, 1000])
ax6.set_yticks([-1000, -500, 0, 500, 1000])

for ax in axs
    ax.tick_params(axis="both", length=ticklength, width=linewidth, 
    color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
end

titles = [L"G_x", L"G_y", "k-space trajectory"]
xlabels = ["Time [ms]", "Time [ms]", L"k_x \ [m^{-1}]"]
ylabels = [L"G_x \ [mT/m]", L"G_y \ [mT/m]", L"k_y \ [m^{-1}]"]
for row = 1:2
    for col = 1:3
        ax = axs[row,col]
        title = titles[col]
        xlabel = xlabels[col]
        ylabel = ylabels[col]
        # ax.set_title(title, fontsize=fontsize_label, color=color_label, pad=6)
        ax.set_xlabel(xlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.set_ylabel(ylabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    end
end

for ax in axs[1, 1:2] 
    ax.set_xlim(t_1p0[1]-1, t_1p0[end]+1)
    ax.set_ylim(-vmax_g_1p0, vmax_g_1p0)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax_g_1p0/4, sigdigits=1)))
end
for ax in axs[2, 1:2] 
    ax.set_xlim(t_0p5[1]-1, t_0p5[end]+1)
    ax.set_ylim(-vmax_g_0p5, vmax_g_0p5)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax_g_0p5/4, sigdigits=1)))
end

# 1p0
for seg = 1:nSegment_1p0
    seg_r = Segment_range_1p0[seg]
    axs[1,1].plot(t_1p0[seg_r], bfieldStitched_1p0[seg_r, 2], color=color_segment[seg], linewidth=linewidth_plot, label="Segment $(seg)")
    axs[1,2].plot(t_1p0[seg_r], bfieldStitched_1p0[seg_r, 3], color=color_segment[seg], linewidth=linewidth_plot, label="Segment $(seg)")
end
axs[1,1].plot(t_1p0, (bfieldStitched_1p0[:, 2] - bfieldStandard_1p0[:, 2]), color=color_difference, linewidth=linewidth_plot, label="Difference")
axs[1,2].plot(t_1p0, (bfieldStitched_1p0[:, 3] - bfieldStandard_1p0[:, 3]), color=color_difference, linewidth=linewidth_plot, label="Difference")
axs[1,3].plot(ksphaStitched_1p0[:,2]/2π, ksphaStitched_1p0[:,3]/2π, color=color_gx, linewidth=linewidth_plot, label="Stitched")
axs[1,3].plot(ksphaStandard_1p0[:,2]/2π, ksphaStandard_1p0[:,3]/2π, color=color_gy, linewidth=linewidth_plot, label="Standard")
axs[1,3].set_aspect(1)

# 0p5
for seg = 1:nSegment_0p5
    nColor = length(color_segment)
    color_idx = seg%nColor != 0 ? seg%nColor : nColor    
    seg_r = Segment_range_0p5[seg]
    axs[2,1].plot(t_0p5[seg_r], bfieldStitched_0p5[seg_r, 2], color=color_segment[color_idx], linewidth=linewidth_plot, label="Stitched")
    axs[2,2].plot(t_0p5[seg_r], bfieldStitched_0p5[seg_r, 3], color=color_segment[color_idx], linewidth=linewidth_plot, label="Stitched")
end 
lDiff_gx, = axs[2,1].plot(t_0p5, (bfieldStitched_0p5[:, 2] - bfieldStandard_0p5[:, 2]), color=color_difference, linewidth=linewidth_plot, label="Difference")
lDiff_gy, = axs[2,2].plot(t_0p5, (bfieldStitched_0p5[:, 3] - bfieldStandard_0p5[:, 3]), color=color_difference, linewidth=linewidth_plot, label="Difference")
axs[2,3].plot(ksphaStitched_0p5[:,2]/2π, ksphaStitched_0p5[:,3]/2π, color=color_gx, linewidth=linewidth_plot, label="Stitched")
axs[2,3].plot(ksphaStandard_0p5[:,2]/2π, ksphaStandard_0p5[:,3]/2π, color=color_gy, linewidth=linewidth_plot, label="Standard")
axs[2,3].set_aspect(1)


# legend
# for (ax, ncol) in zip(axs[[5,6]], [2,2])   # setup legend for each subplot
#     ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
#     loc="center left", bbox_to_anchor=(0,1.0),
#     frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
# end
for (ax, color_l) in zip(axs[[5,6]], color_segment[1:2])
    ax.legend(
        [lines.Line2D([], [], color=color_segment[1], linewidth=linewidth_legend), 
            lines.Line2D([], [], color=color_segment[2], linewidth=linewidth_legend)], 
        ["Stitched", "Standard"], 
        handler_map=Dict(Tuple(ls) => legend_handler.HandlerTuple(ndivide=length(color_segment), pad=0)),
        fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
        loc="center left", bbox_to_anchor=(0,1.0),
        frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

for (ax, ncol) in zip(axs[[1,3]], [2,2])   # setup legend for each subplot
    ls = []
    for color in color_segment[1:4]
        push!(ls, lines.Line2D([], [], color=color, linewidth=linewidth_legend))
    end
    ax.legend([Tuple(ls), lines.Line2D([], [], color=color_difference, linewidth=linewidth_legend)], ["Stitched", "Difference"], 
        handler_map=Dict(Tuple(ls) => legend_handler.HandlerTuple(ndivide=4, pad=0)),
        fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
        loc="center left", bbox_to_anchor=(0,1.0),
        frameon=false, handlelength=2.0, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

ls = []
for color in color_segment
    push!(ls, lines.Line2D([], [], color=color, linewidth=linewidth_legend))
end
for (ax, lDiff) in zip(axs[2,1:2], [lDiff_gx, lDiff_gy])
    ax.legend([Tuple(ls), lines.Line2D([], [], color=color_difference, linewidth=linewidth_legend)], ["Stitched", "Difference"], 
        handler_map=Dict(Tuple(ls) => legend_handler.HandlerTuple(ndivide=length(color_segment), pad=0)),
        fontsize=fontsize_legend, labelcolor=color_label, ncols=2, 
        loc="center left", bbox_to_anchor=(0,1.0),
        frameon=false, handlelength=2.0, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end


fig.align_ylabels()
orders = ["a", "b", "c", "d", "e", "f"]
for row = 1:2
    for col = 1:3
        order = orders[3*(row-1)+col]
        fig.text(0.01+(col-1)/3, 1-0.5*(row-1), "($(order))", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)
    end
end

fig.tight_layout(pad=0, h_pad=0.5, w_pad=0.1)
fig.savefig("$(outpath)/Fig7.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig7.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(@__DIR__)/Figures/tiff/Fig7.tiff", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
