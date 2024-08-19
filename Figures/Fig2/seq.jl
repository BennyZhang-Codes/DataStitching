using PyPlot
using KomaHighOrder

sh = SphericalHarmonics()
nSegment_fully = 4
nSegment_under = 36

_, ntStitched_fully, _ = load_dfc(;dfc_method=:Stitched, seqname="demo", r=1);
_, ntStitched_under, _ = load_dfc(;dfc_method=:Stitched, seqname="demo", r=30);

hoseq_fully = demo_hoseq(dfc_method=:Stitched)[4:end]   
hoseq_under = demo_hoseq(dfc_method=:Stitched, r=30)[4:end]   


samples_fully = get_samples(hoseq_fully; off_val=Inf);
samples_under = get_samples(hoseq_under; off_val=Inf);

t_adc_fully = KomaMRIBase.get_adc_sampling_times(hoseq_fully.SEQ);
t_adc_under = KomaMRIBase.get_adc_sampling_times(hoseq_under.SEQ);

dt = 1e-6
a = Int64.(vec(ntStitched_fully)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger_fully = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq_fully.SEQ)[5]

a = Int64.(vec(ntStitched_under)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger_under = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq_under.SEQ)[5]


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 18/2.54
figure_height      = 5/2.54
linewidth          = 0.5
ticklength         = 1.5
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
fontsize_subfigure = 8
pad_labeltick      = 2
pad_label          = 2
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_trigger      = "#000000"

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(figure_width, figure_height), facecolor=color_facecoler)

samples_list = [samples_fully, samples_under]
t_trigger_list = [t_trigger_fully, t_trigger_under]
t_adc_list = [t_adc_fully, t_adc_under]



for (ax, samples, t_trigger) in zip(axs, samples_list, t_trigger_list)
    ax.tick_params(axis="both", length=ticklength, width=linewidth, color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
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

    ax.scatter(t_trigger*1e3, -60*ones(length(t_trigger)), s=15, marker=L"\mapsup", color=color_trigger, linewidth=0.1, label=L"trigger")

    ax.set_ylim(-75, 75)
    ax.set_ylabel("amplitude (mT/m, mG)", fontsize=fontsize_label, color=color_label)
    ax.set_xlabel("time (ms)", fontsize=fontsize_label, color=color_label)
    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=5, loc="upper left", frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)
end

fig.align_ylabels()
orders = ["a", "b"]

for col = 1:2
    order = orders[col]
    fig.text(0.01+(col-1)/2, 1, "($(order))", ha="left", va="baseline", fontsize=fontsize_subfigure)
end
fig.tight_layout(pad=0, h_pad=0.5, w_pad=0.5)
fig.savefig("Figures/Fig2/Fig2.png", dpi=300, bbox_inches="tight")