using PyPlot
using KomaHighOrder

key = :Standard
hoseq_Standard = demo_hoseq(dfc_method=:Standard)[8]   # :Standard or :Stitched
hoseq_Stitched = demo_hoseq(dfc_method=:Stitched)[8]   # :Standard or :Stitched
sh = SphericalHarmonics()

index = [0, 3, 4, 5, 6, 7, 8]

samples_Standard = get_samples(hoseq_Standard; off_val=0) 
samples_Stitched = get_samples(hoseq_Stitched; off_val=0) 
t_adc = samples_Stitched.h0.t .* 1e3  # convert to ms  


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 18/2.54
figure_height      = 9/2.54
linewidth          = 0.3
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_segment      = ["C0", "C1", "C2", "C3"]

fig, axs = plt.subplots(nrows=7, ncols=2, figsize=(figure_width, figure_height), 
                        sharex=true, facecolor=color_facecoler)

for ax in axs
    ax.tick_params(axis="both", length=linewidth*5, width=linewidth, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecoler)
end
ylabels = [L"G_{0}", L"G_{x}", L"G_{y}", L"G_{z}", L"G_{xy}", L"G_{zy}", L"G_{2z^2-(x^2+y^2)}", L"G_{xz}", L"G_{x^2-y^2}"]

vmax = [0.1, 2, 5, 5, 5, 5, 5]
for idx in eachindex(index)
    ax1, ax2 = axs[idx, :]
    term = index[idx]
    vlim = vmax[idx]
    println(vlim)
    ax1.plot(t_adc, (samples_Stitched[Symbol("h$(term)")].A-samples_Standard[Symbol("h$(term)")].A)*1e3, color="black", linewidth=linewidth, label="Difference")
    ax1.plot(t_adc, samples_Stitched[Symbol("h$(term)")].A*1e3, color=sh.dict["h$(term)"].color, linewidth=linewidth, label="Stitched measurement")
    ax1.set_ylim(-vlim, vlim)
    ax1.set_ylabel("$(ylabels[term+1])\n(m$(sh.dict["h$(term)"].unit))",rotation=0, 
        ha="right", va="center", x=0, y=0.5,
        fontsize=fontsize_label, color=color_label)

    ax2.plot(t_adc, (samples_Stitched[Symbol("h$(term)")].A-samples_Standard[Symbol("h$(term)")].A)*1e3, color="black", linewidth=linewidth, label="Difference")
    ax2.plot(t_adc, samples_Stitched[Symbol("h$(term)")].A*1e3, color=sh.dict["h$(term)"].color, linewidth=linewidth, label="Stitched measurement")
    ax2.set_ylim(-vlim, vlim)
    ax2.set_ylabel("$(ylabels[term+1])\n(m$(sh.dict["h$(term)"].unit))",rotation=0, 
        ha="right", va="center", x=0, y=0.5,
        fontsize=fontsize_label, color=color_label)
end

axs[7, 1].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label)
axs[7, 2].set_xlabel("Time (ms)", fontsize=fontsize_label, color=color_label)

# fig.align_ylabels()
for (ax, ncol) in zip(axs, ones(14)*2)   # setup legend for each subplot
    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
    loc="upper left", 
    bbox_to_anchor=(0,1.05),
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end

fig.align_ylabels()
fig.tight_layout(pad=0, h_pad=0.1)
fig.savefig("Figures/Fig4/Fig4.png", dpi=300, bbox_inches="tight")

