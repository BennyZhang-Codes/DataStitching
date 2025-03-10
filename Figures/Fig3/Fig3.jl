using PyPlot
using HighOrderMRI

outpath = "$(@__DIR__)/Figures/Fig3/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory

sh = SphericalHarmonics()
nSegment_fully = 4
nSegment_under = 36


seq_img_file_1p0 = "$(@__DIR__)/Figures/Fig3/7T_1p0_200_r4_max51_fa90.seq" 
seq_dfc_file_1p0 = "$(@__DIR__)/Figures/Fig3/7T_1p0_200_r4_seg4_max51.seq" 

seq_img_file_0p5 = "$(@__DIR__)/Figures/Fig3/7T_0p5_400_r4_max51_fa90.seq" 
seq_dfc_file_0p5 = "$(@__DIR__)/Figures/Fig3/7T_0p5_400_r4_seg276_max51.seq" 

seq_1p0 = read_seq(seq_img_file_1p0)[2:8]
seq_0p5 = read_seq(seq_img_file_0p5)[2:8]
seq_1p0.GR[1,:] = -seq_1p0.GR[1,:]
seq_0p5.GR[1,:] = -seq_0p5.GR[1,:]


trigger_delays_1p0 = read_seq(seq_dfc_file_1p0).DEF["skope_triggerDelays"]
trigger_delays_0p5 = read_seq(seq_dfc_file_0p5).DEF["skope_triggerDelays"]

samples_1p0 = get_samples(seq_1p0; off_val=Inf);
samples_0p5 = get_samples(seq_0p5; off_val=Inf);

t_adc_1p0 = KomaMRIBase.get_adc_sampling_times(seq_1p0);
t_adc_0p5 = KomaMRIBase.get_adc_sampling_times(seq_0p5);

dt = 1e-6
delay = 0.0005
# a = Int64.(vec(ntStitched_fully)); a[1] += 3; e = cumsum(a); s = e .- a;
# t_trigger_1p0 = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(seq_1p0)[5] .- delay

# a = Int64.(vec(ntStitched_under)); a[1] += 3; e = cumsum(a); s = e .- a;
# t_trigger_0p5 = (s .+ 1)*dt .+ KomaMRIBase.get_block_start_times(seq_0p5)[5] .- delay


t_trigger_1p0 = trigger_delays_1p0 .+ KomaMRIBase.get_block_start_times(seq_1p0)[6] .- delay
t_trigger_0p5 = trigger_delays_0p5 .+ KomaMRIBase.get_block_start_times(seq_0p5)[6] .- delay




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



# fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw=Dict("width_ratios"=> [samples_1p0.gx.t[end], samples_0p5.gx.t[end]]), figsize=(figure_width, figure_height), facecolor=color_facecolor)
fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw=Dict("width_ratios"=> [4,6]), figsize=(figure_width, figure_height), facecolor=color_facecolor)

samples_list = [samples_1p0, samples_0p5]
t_trigger_list = [t_trigger_1p0, t_trigger_0p5]
t_adc_list = [t_adc_1p0, t_adc_0p5]

axs[1].set_xticks([0,10,20,30,40])
axs[2].set_xticks([0,10,20,30,40,50,60,70,80,90,100])

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

    ax.plot(samples.rf.t*1e3, abs.(samples.rf.A)*1e6, color=color_rf, linewidth=linewidth, label=L"|B_{1}|")

    # ax.scatter(t_trigger*1e3, -60*ones(length(t_trigger)), s=markersize_trigger, marker=L"\uparrow", color=color_trigger, linewidth=linewidth_marker, label=L"Trigger")
    ax.scatter(t_trigger*1e3, -60*ones(length(t_trigger)), s=markersize_trigger, marker="|", color=color_trigger, linewidth=linewidth_marker, label=L"Trigger")

    ax.set_ylim(-75, 75)
    ax.set_ylabel("Amplitude [mT/m, uT]", fontsize=fontsize_label, color=color_label)
    ax.set_xlabel("Time [ms]", fontsize=fontsize_label, color=color_label)
    ax.legend(loc="center left", bbox_to_anchor=(0, 1.02), fontsize=fontsize_legend, labelcolor=color_label, 
        scatteryoffsets=[0.5], borderpad=0,
        ncols=5, frameon=false, handlelength=1, handletextpad=0.2, columnspacing=0.5)

    # for trigger in t_trigger
    #     ax.arrow(trigger*1e3, -65, 0, 10, overhang=0.1, width=0.1, head_width=0.2, head_length=3, linewidth=0.1, length_includes_head=true, fc=color_label)
    # end
end

# fig.legend(loc="center", bbox_to_anchor=(0.5, 1.1),
#     handles=[axs[1].lines[2], axs[1].lines[3], axs[1].lines[4], axs[1].lines[5], legend_trigger], 
#     fontsize=fontsize_legend, labelcolor=color_label, 
#     scatteryoffsets=[0.5],
#     ncols=5, frameon=false, handlelength=1, handletextpad=0.2, columnspacing=0.5)

fig.align_ylabels()

fig.text(0, 1, "(a)", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)
# fig.text(0.349, 1, "(b)", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)
fig.text(0.415, 1, "(b)", ha="left", va="center", fontsize=fontsize_subfigure, color=color_label)

fig.tight_layout(pad=0, h_pad=0, w_pad=0.5)

fig.savefig("$(outpath)/Fig3.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(outpath)/Fig3.svg", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
fig.savefig("$(@__DIR__)/Figures/tiff/Fig3.tiff", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)




