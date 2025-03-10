using PyPlot
using HighOrderMRI
lines = PyPlot.matplotlib.lines

outpath = "$(@__DIR__)/Figures/Fig5/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory


sh = SphericalHarmonics()
terms = [0, 3, 4, 5, 6, 7, 8]

hoseq_fully_Standard = load_hoseq(dfc_method=:Standard)[8];   # :Standard or :Stitched
hoseq_fully_Stitched = load_hoseq(dfc_method=:Stitched)[8];   # :Standard or :Stitched
hoseq_under_Standard = load_hoseq(dfc_method=:Standard, r=30)[8];   # :Standard or :Stitched
hoseq_under_Stitched = load_hoseq(dfc_method=:Stitched, r=30)[8];   # :Standard or :Stitched


samples_fully_Standard = get_samples(hoseq_fully_Standard; off_val=0); 
samples_fully_Stitched = get_samples(hoseq_fully_Stitched; off_val=0); 
t_adc_fully = samples_fully_Stitched.h0.t .* 1e3;  # convert to ms  
samples_under_Standard = get_samples(hoseq_under_Standard; off_val=0); 
samples_under_Stitched = get_samples(hoseq_under_Stitched; off_val=0); 
t_adc_under = samples_under_Stitched.h0.t .* 1e3;  # convert to ms  



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

fig, axs = plt.subplots(nrows=7, ncols=2, figsize=(figure_width, figure_height), 
                        facecolor=color_facecolor)

for ax in axs[1:6, :]
    ax.tick_params(axis="both", labelbottom=false)
end
axs[7,1].set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80])
axs[7,2].set_xticks([0, 5, 10, 15, 20])
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


vmaxs = [0.15 0.15; 1.5 1.5; 8 5; 8 5; 3 3; 10 6; 4 4];
t_adc = [t_adc_fully, t_adc_under];
samples_Stitched = [samples_fully_Stitched, samples_under_Stitched];
samples_Standard = [samples_fully_Standard, samples_under_Standard];
ncols = ones(14)*2;
for row = 1:7
    for col = 1:2
        ax = axs[row, col]
        term = terms[row]
        vmax = vmaxs[row, col]
        t = t_adc[col]
        sampleStitched = samples_Stitched[col]
        sampleStandard = samples_Standard[col]
        ncol = ncols[row]
        ax.set_ylim(-vmax, vmax)
        ax.set_xlim(t[1], t[end])

        line1, = ax.plot(t, (sampleStitched[Symbol("h$(term)")].A-sampleStandard[Symbol("h$(term)")].A)*1e3, color=color_difference, linewidth=linewidth_plot, label="Difference")
        line2, = ax.plot(t, sampleStitched[Symbol("h$(term)")].A*1e3, color=sh.dict["h$(term)"].color, linewidth=linewidth_plot, label="Stitched")
        ax.yaxis.set_major_locator(plt.MultipleLocator(vmax))
        ax.set_ylabel("$(ylabels[term+1])\n[$(units[term+1])]",rotation=0, 
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

fig.savefig("$(outpath)/Fig5.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig5.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(@__DIR__)/Figures/tiff/Fig5.tiff", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)

