using PyPlot
using HighOrderMRI

outpath = "$(@__DIR__)/Figures/Fig2/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

sh = SphericalHarmonics()
nSegment_1p0 = 4
nSegment_0p3 = 36

_, ntStitched_1p0, _ = load_dfc(;dfc_method=:Stitched, seqname="spiral", r=1);
_, ntStitched_0p3, _ = load_dfc(;dfc_method=:Stitched, seqname="spiral", r=30);

hoseq_1p0 = load_hoseq(dfc_method=:Stitched)[4:end]   
hoseq_0p3 = load_hoseq(dfc_method=:Stitched, r=30)[4:end]   


samples_1p0 = get_samples(hoseq_1p0; off_val=Inf);
samples_0p3 = get_samples(hoseq_0p3; off_val=Inf);

t_adc_1p0 = KomaMRIBase.get_adc_sampling_times(hoseq_1p0.SEQ);
t_adc_0p3 = KomaMRIBase.get_adc_sampling_times(hoseq_0p3.SEQ);

dt = 1e-6
delay = 0.0005
a = Int64.(vec(ntStitched_1p0)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger_1p0 = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq_1p0.SEQ)[5] .- delay

a = Int64.(vec(ntStitched_0p3)); a[1] += 3; e = cumsum(a); s = e .- a;
t_trigger_0p3 = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(hoseq_0p3.SEQ)[5] .- delay


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 17/2.53999863
figure_height      = 4/2.53999863
linewidth          = 0.6
linewidth_marker   = 0.35
ticklength         = 1.5

fontsize_subfigure = 9
fontsize_label     = 7
fontsize_legend    = 7
fontsize_ticklabel = 6

pad_labeltick      = 2
pad_label          = 2
markersize_trigger = 30

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

# fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw=Dict("width_ratios"=> [samples_1p0.gx.t[end], samples_0p3.gx.t[end]]), figsize=(figure_width, figure_height), facecolor=color_facecolor)
fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw=Dict("width_ratios"=> [6,4]), figsize=(figure_width, figure_height), facecolor=color_facecolor)


samples_list = [samples_1p0, samples_0p3]
t_trigger_list = [t_trigger_1p0, t_trigger_0p3]
t_adc_list = [t_adc_1p0, t_adc_0p3]

axs[1].set_xticks([0,10,20,30,40,50,60,70,80,90])
axs[2].set_xticks([0,10,20,30])

legend_trigger = nothing
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

    # legend_trigger = ax.scatter(t_trigger*1e3, -60*ones(length(t_trigger)), s=markersize_trigger, marker=L"\uparrow", color=color_trigger, linewidth=linewidth_marker, label=L"Trigger")
    legend_trigger = ax.scatter(t_trigger*1e3, -60*ones(length(t_trigger)), s=markersize_trigger, marker="|", color=color_trigger, linewidth=linewidth_marker, label=L"Trigger")

    ax.set_ylim(-75, 75)
    ax.set_ylabel("Amplitude [mT/m, mG]", fontsize=fontsize_label, color=color_label)
    ax.set_xlabel("Time [ms]", fontsize=fontsize_label, color=color_label)
end

# fig.legend(loc="center", bbox_to_anchor=(0.5, 1.05),
#     handles=[axs[1].lines[2], axs[1].lines[3], axs[1].lines[4], axs[1].lines[5], legend_trigger], 
#     fontsize=fontsize_legend, labelcolor=color_label, 
#     scatteryoffsets=[0.5],
#     ncols=5, frameon=false, handlelength=1, handletextpad=0.2, columnspacing=0.5)

axs[1].legend(loc="center left", bbox_to_anchor=(0, 1.02), fontsize=fontsize_legend, labelcolor=color_label, 
scatteryoffsets=[0.5], borderpad=0,
ncols=5, frameon=false, handlelength=1, handletextpad=0.2, columnspacing=0.5)
axs[2].legend(loc="center left", bbox_to_anchor=(0, 1.02), fontsize=fontsize_legend, labelcolor=color_label, 
scatteryoffsets=[0.5], borderpad=0,
ncols=5, frameon=false, handlelength=1, handletextpad=0.2, columnspacing=0.5)

fig.align_ylabels()
orders = ["a", "b"]

fig.text(0, 1, "(a)", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)
# fig.text(0.725, 1, "(b)", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)
fig.text(0.595, 1, "(b)", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)

fig.tight_layout(pad=0, h_pad=0, w_pad=0.5)

fig.savefig("$(outpath)/Fig2.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig2.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(@__DIR__)/Figures/tiff/Fig2.tiff", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)


