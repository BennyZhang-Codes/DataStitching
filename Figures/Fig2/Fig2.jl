using PyPlot
using KomaHighOrder

sh = SphericalHarmonics()
nSegment_fully = 4
nSegment_under = 36

_, ntStitched_fully, _ = load_dfc(;dfc_method=:Stitched, seqname="spiral", r=1);
_, ntStitched_under, _ = load_dfc(;dfc_method=:Stitched, seqname="spiral", r=30);

hoseq_fully = load_hoseq(dfc_method=:Stitched)[4:end]   
hoseq_under = load_hoseq(dfc_method=:Stitched, r=30)[4:end]   


samples_fully = get_samples(hoseq_fully; off_val=Inf);
samples_under = get_samples(hoseq_under; off_val=Inf);

t_adc_fully = KomaMRIBase.get_adc_sampling_times(hoseq_fully.SEQ);
t_adc_under = KomaMRIBase.get_adc_sampling_times(hoseq_under.SEQ);

dt = 1e-6
delay = 0.0005
a = Int64.(vec(ntStitched_fully)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger_fully = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq_fully.SEQ)[5] .- delay

a = Int64.(vec(ntStitched_under)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger_under = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq_under.SEQ)[5] .- delay


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 17/2.53999863
figure_height      = 4/2.53999863
linewidth          = 0.8
linewidth_marker   = 0.5
ticklength         = 1.5
fontsize_legend    = 7
fontsize_label     = 7
fontsize_ticklabel = 6
fontsize_subfigure = 9
pad_labeltick      = 2
pad_label          = 2
markersize_trigger = 40

color_facecolor    = "#ffffff"
color_label        = "#000000"
color_trigger      = "#555555"

color_gx           = "#636EFA"
color_gy           = "#EF553B"
color_gz           = "#00CC96"
color_rf           = "#AB63FA"

color_gx           = "#6699CC"
color_gy           = "#FF6666"
color_gz           = "#669966"
color_rf           = "#CC66CC"

# color_gx           = "C0"
# color_gy           = "C3"
# color_gz           = "C2"
# color_rf           = "C5"

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(figure_width, figure_height), facecolor=color_facecolor)

samples_list = [samples_fully, samples_under]
t_trigger_list = [t_trigger_fully, t_trigger_under]
t_adc_list = [t_adc_fully, t_adc_under]


for (ax, samples, t_trigger) in zip(axs, samples_list, t_trigger_list)
    ax.tick_params(axis="both", length=ticklength, width=linewidth, color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_xlim(samples.gz.t[1]*1e3 - 1, samples.gz.t[end]*1e3+1)
    ax.set_facecolor(color_facecolor)
    ax.yaxis.set_major_locator(plt.MultipleLocator(20))

    ax.axhline(y=0, color="#555555", linewidth=linewidth)
    ax.plot(samples.gx.t*1e3,       samples.gx.A*1e3, color=color_gx, linewidth=linewidth, label=L"G_{x}")
    ax.plot(samples.gy.t*1e3,       samples.gy.A*1e3, color=color_gy, linewidth=linewidth, label=L"G_{y}")
    ax.plot(samples.gz.t*1e3,       samples.gz.A*1e3, color=color_gz, linewidth=linewidth, label=L"G_{z}")

    ax.plot(samples.rf.t*1e3, abs.(samples.rf.A)*1e7, color=color_rf, linewidth=linewidth, label=L"|B_{1}|")

    ax.scatter(t_trigger*1e3, -60*ones(length(t_trigger)), s=markersize_trigger, marker=L"\uparrow", color=color_trigger, linewidth=linewidth_marker, label=L"Trigger")

    ax.set_ylim(-75, 75)
    ax.set_ylabel("Amplitude [mT/m, mG]", fontsize=fontsize_label, color=color_label)
    ax.set_xlabel("Time [ms]", fontsize=fontsize_label, color=color_label)
    ax.legend(bbox_to_anchor=(0.01, 1.1), fontsize=fontsize_legend, labelcolor=color_label, 
        scatteryoffsets=[0.5],
        ncols=5, loc="upper left", frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)
end

fig.align_ylabels()
orders = ["a", "b"]

for col = 1:2
    order = orders[col]
    fig.text(0+(col-1)/2, 1, "($(order))", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)
end
fig.tight_layout(pad=0, h_pad=0.5, w_pad=0.8)
fig.savefig("Figures/Fig2/Fig2.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("Figures/Fig2/Fig2.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)