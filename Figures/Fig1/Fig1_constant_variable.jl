using PyPlot
Line2D = matplotlib.lines.Line2D
Rectangle = matplotlib.patches.Rectangle
Patch = matplotlib.patches.Patch

outpath = "$(@__DIR__)/Figures/Fig1/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 15/2.53999863
figure_height      = 8/2.53999863

fontsize_legend    = 8
fontsize_title     = 9
fontsize_label     = 7
markersize_trigger = 80

color_facecolor      = "#ffffff"
color_title          = "#000000"
color_label          = "#555555"
color_text           = "#ffffff"
color_G              = "#6699CC"
color_F              = "#FF6666"
color_seg            = "#99CC99"
color_trigger        = "#555555"

linewidth            = 0.8
linewidth_marker     = 0.5

fig, axs = plt.subplots(2, 1, figsize=(figure_width, figure_height), facecolor=color_facecolor)
y_G = 0.5
y_F = 0.3
h_G = 0.15
h_F = 0.1


legend_trigger = nothing
modes = ["Constant-segment stitching", "Variable-segment stitching"]
for (ax, mode, duration_segs) in zip(axs, modes, [[15, 15, 15, 15], [25, 20, 10, 5]])
    ax.set_title(mode, fontsize=fontsize_title, color=color_title)
    # sequence 
    nTR = 4
    TR = 90
    TE = 5
    duration = 60
    TRStarts = collect(range(0,nTR*TR, nTR+1))[1:nTR]
    ReadoutStarts = TRStarts .+ TE
    # segments
    # duration_segs = [25, 20, 10, 5]
    SegEnds = cumsum(duration_segs)
    SegStarts = SegEnds .- duration_segs
    # field camera
    delay = 2
    CameraStarts = SegStarts .- delay
    CameraEnds = CameraStarts .+ maximum(duration_segs) .+ delay .+ 1

    # 标签
    ax.text(-8, y_G, "G", fontsize=fontsize_label, color=color_G, verticalalignment="center", horizontalalignment="right")
    ax.text(-8, y_F, "F", fontsize=fontsize_label, color=color_F, verticalalignment="center", horizontalalignment="right")
    # 绘制水平线 G 和 F
    ax.add_line(Line2D([-5, nTR*TR+10], [y_F, y_F], linewidth=linewidth, color=color_F))
    ax.add_line(Line2D([-5, nTR*TR+10], [y_G, y_G], linewidth=linewidth, color=color_G))

    ax.set_ylim(y_F-0.2, y_G+0.35)
    ax.set_xlim(-5, nTR*TR+15)

    ax.text(nTR*TR+10, y_G+0.1, "...", fontsize=fontsize_label, color=color_label, verticalalignment="center", horizontalalignment="right")
    ax.text(nTR*TR+10, y_F+0.1, "...", fontsize=fontsize_label, color=color_label, verticalalignment="center", horizontalalignment="right")

    ax.axis("off")

    for idx_tr = 1:nTR
        seg_s    = ReadoutStarts[idx_tr] + SegStarts[idx_tr]
        seg_e    = ReadoutStarts[idx_tr] + SegEnds[idx_tr]
        camera_s = ReadoutStarts[idx_tr] + CameraStarts[idx_tr]
        camera_e = ReadoutStarts[idx_tr] + CameraEnds[idx_tr]
        ax.add_patch(Rectangle([camera_s,    y_F], camera_e-camera_s,  h_F, facecolor=color_F,   linewidth=linewidth, edgecolor=color_F))
        ax.add_patch(Rectangle([   seg_s,    y_F],       seg_e-seg_s,  h_F, facecolor=color_seg, linewidth=linewidth, edgecolor=color_seg, hatch="/////", fill=false, zorder=99))

    end

    for idx_tr = 1:nTR
        ax.add_line(Line2D([TRStarts[idx_tr], TRStarts[idx_tr]], [y_F, y_G+0.35], linewidth=linewidth, color=color_label, linestyle="--"))
        if idx_tr == nTR
            ax.add_line(Line2D([TRStarts[idx_tr]+TR, TRStarts[idx_tr]+TR], [y_F, y_G+0.35], linewidth=linewidth, color=color_label, linestyle="--"))
        end
        pos = ReadoutStarts[idx_tr]
        # ax.add_line(Line2D([         pos, pos+duration], [1.3, 1.3], linewidth=2, color="#000000"))
        # ax.add_line(Line2D([         pos,          pos], [  1, 1.3], linewidth=2, color="#000000"))
        # ax.add_line(Line2D([pos+duration, pos+duration], [  1, 1.3], linewidth=2, color="#000000"))
        ax.fill_between([pos,pos+duration],  y_G, y_G+h_G, color=color_G, alpha=1, linewidth=linewidth)

        ax.text(pos+duration/2, (y_G+y_G+h_G)/2,                     "Readout", fontsize=fontsize_label, color=color_text,  horizontalalignment="center", verticalalignment="center")
        ax.text(TRStarts[idx_tr]+TR/2,  y_G+h_G+0.1, "TR #$(idx_tr)", fontsize=fontsize_label, color=color_label, horizontalalignment="center", verticalalignment="center")
        ax.arrow(   duration/3+TRStarts[idx_tr], y_G+h_G+0.1, -duration/3, 0, overhang=0.3, width=0.005, head_width=0.025, head_length=5, linewidth=0, length_includes_head=true, fc=color_label)
        ax.arrow(TR-duration/3+TRStarts[idx_tr], y_G+h_G+0.1,  duration/3, 0, overhang=0.3, width=0.005, head_width=0.025, head_length=5, linewidth=0, length_includes_head=true, fc=color_label)
        for idx = 1:nTR
            seg_s = pos + SegStarts[idx]
            seg_e = pos + SegEnds[idx]
            if idx == idx_tr
                ax.add_line(Line2D([seg_s, seg_s], [y_F+h_F, y_G], linewidth=linewidth, color=color_G, linestyle="--"))
                ax.add_line(Line2D([seg_e, seg_e], [y_F+h_F, y_G], linewidth=linewidth, color=color_G, linestyle="--"))
                # ax.arrow(seg_s, y_G-0.05, 0, 0.05, overhang=0, width=0.5, head_width=2, head_length=0.05, linewidth=0, length_includes_head=true, fc="C2")
                # ax.arrow(seg_e, y_G-0.05, 0, 0.05, overhang=0, width=0.5, head_width=2, head_length=0.05, linewidth=0, length_includes_head=true, fc="C2")
            end
            # if idx < idx_tr
            #     ax.add_line(Line2D([seg_s, seg_s], [y_G-0.03, y_G], linewidth=linewidth, color=color_G))
            # end
            # if idx > idx_tr
            #     ax.add_line(Line2D([seg_e, seg_e], [y_G-0.03, y_G], linewidth=linewidth, color=color_G))
            # end
        end
    end
    legend_trigger = ax.scatter(CameraStarts+ReadoutStarts, ones(nTR)*0.25, s=markersize_trigger, marker=L"\uparrow", color=color_trigger, linewidth=linewidth_marker, label="Trigger")
end

fig.legend(handles=[
        legend_trigger,
        Patch(color=color_F, label="Recording window"),
        # Patch(color=color_seg, label="Effective data segment"),
        Rectangle([0,0],0,  0, label="Effective data segment",facecolor=color_seg, linewidth=linewidth, edgecolor=color_seg, hatch="/////", fill=false)], 
    ncol=3, bbox_to_anchor=(0.5, 0.0), frameon=false, fontsize=fontsize_legend, loc="lower center", 
    scatteryoffsets=[0.7],
    columnspacing=1, borderpad=0, handletextpad=0.3, handleheight=0.7, handlelength=1.5)
fig.tight_layout(h_pad=0)

fig.savefig("$(outpath)/Fig1.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig1.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
