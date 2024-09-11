using PyPlot
using KomaHighOrder

sh = SphericalHarmonics()
terms = [0, 3, 4, 5, 6, 7, 8]

hoseq_fully_Standard = demo_hoseq(dfc_method=:Standard)[8];   # :Standard or :Stitched
hoseq_fully_Stitched = demo_hoseq(dfc_method=:Stitched)[8];   # :Standard or :Stitched
hoseq_under_Standard = demo_hoseq(dfc_method=:Standard, r=30)[8];   # :Standard or :Stitched
hoseq_under_Stitched = demo_hoseq(dfc_method=:Stitched, r=30)[8];   # :Standard or :Stitched


samples_fully_Standard = get_samples(hoseq_fully_Standard; off_val=0); 
samples_fully_Stitched = get_samples(hoseq_fully_Stitched; off_val=0); 
t_adc_fully = samples_fully_Stitched.h0.t .* 1e3;  # convert to ms  
samples_under_Standard = get_samples(hoseq_under_Standard; off_val=0); 
samples_under_Stitched = get_samples(hoseq_under_Stitched; off_val=0); 
t_adc_under = samples_under_Stitched.h0.t .* 1e3;  # convert to ms  



# matplotlib.rc("mathtext", default="regular")
# matplotlib.rc("figure", dpi=200)
# matplotlib.rc("font", family="Times New Roman")
# matplotlib.rcParams["mathtext.default"]
figure_width       = 18/2.54
figure_height      = 9/2.54
# linewidth          = 0.5
# ticklength         = 1.5
# fontsize_legend    = 5
# fontsize_label     = 6
# fontsize_ticklabel = 4
# fontsize_subfigure = 8
# pad_labeltick      = 2
# pad_label          = 2
# color_facecoler    = "#ffffff"
# color_label        = "#000000"
color_segment      = ["C0", "C1", "C2", "C3"]

fig, axs = plt.subplots(nrows=7, ncols=2, figsize=(figure_width, figure_height), 
                        facecolor=color_facecoler)

for ax in axs[1:6, :]
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
    ax.set_facecolor(color_facecoler)
end
ylabels = [L"G_{0}", L"G_{x}", L"G_{y}", L"G_{z}", L"G_{xy}", L"G_{zy}", L"G_{2z^2-(x^2+y^2)}", L"G_{xz}", L"G_{x^2-y^2}"]

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

        line1, = ax.plot(t, (sampleStitched[Symbol("h$(term)")].A-sampleStandard[Symbol("h$(term)")].A)*1e3, color="black", linewidth=linewidth, label="Difference")
        line2, = ax.plot(t, sampleStitched[Symbol("h$(term)")].A*1e3, color=sh.dict["h$(term)"].color, linewidth=linewidth, label="Stitched measurement")
        ax.yaxis.set_major_locator(plt.MultipleLocator(vmax))
        ax.set_ylabel("$(ylabels[term+1])\n(m$(sh.dict["h$(term)"].unit))",rotation=0, 
            ha="right", va="center", x=0, y=0.5,
            fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.legend(handles=[line2, line1], fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
            loc="upper left", bbox_to_anchor=(0,1.05),
            frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
    end
end


axs[7, 1].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)
axs[7, 2].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label, labelpad=pad_label)

fig.text(0.01, 1, "(a)", ha="left", va="baseline", fontsize=fontsize_subfigure)
fig.text(0.51, 1, "(b)", ha="left", va="baseline", fontsize=fontsize_subfigure)

fig.align_ylabels()
# fig.subplots_adjust(left=0.1, right=0.95, bottom=0.05, top=0.95, wspace=0.2, hspace=0.2)
fig.tight_layout(pad=0, h_pad=0.1, w_pad=0.3)
# fig.savefig("Figures/Fig4/Fig4.png", dpi=300, bbox_inches="tight")

