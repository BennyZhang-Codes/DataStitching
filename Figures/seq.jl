using PyPlot
using KomaHighOrder


hoseq = demo_hoseq(skope_method=key)   # :Standard or :Stitched 
samples = get_samples(hoseq; off_val=Inf)


figure_width       = 9
figure_height      = 5
linewidth          = 1
fontsize_legend    = 10
fontsize_label     = 11
fontsize_ticklabel = 8
color_facecoler    = "#ffffff"
color_label        = "#000000"

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecoler)

ax.tick_params(axis="both", color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_color(color_label)
    ax.spines[spine].set_visible(false)
end
ax.set_xlim(samples.gz.t[1]*1e3, samples.gz.t[end]*1e3)
ax.set_facecolor(color_facecoler)
ax.yaxis.set_major_locator(plt.MultipleLocator(20))

ax.axhline(y=0, color="#555555", linewidth=linewidth)
ax.plot(samples.gx.t*1e3,       samples.gx.A*1e3, color="#636EFA", linewidth=linewidth, label=L"G_{x}")
ax.plot(samples.gy.t*1e3,       samples.gy.A*1e3, color="#EF553B", linewidth=linewidth, label=L"G_{y}")
ax.plot(samples.gz.t*1e3,       samples.gz.A*1e3, color="#00CC96", linewidth=linewidth, label=L"G_{z}")

ax.plot(samples.rf.t*1e3, abs.(samples.rf.A)*1e7, color="#AB63FA", linewidth=linewidth, label=L"|B_{1}|")

ax.set_ylim(-80, 120)

ax.set_ylabel("amplitude (mT/m, mG)", fontsize=fontsize_label, color=color_label)
ax.set_xlabel("time (ms)", fontsize=fontsize_label, color=color_label)
ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=4, loc="upper right", frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)
fig.tight_layout()

fig.savefig("Figures/out/seq.png", dpi=300)