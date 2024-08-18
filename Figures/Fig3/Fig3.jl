using PyPlot
using KomaHighOrder

key = :Standard
hoseq_Standard = demo_hoseq(dfc_method=:Standard)[8]   # :Standard or :Stitched
hoseq_Stitched = demo_hoseq(dfc_method=:Stitched)[8]   # :Standard or :Stitched
sh = SphericalHarmonics()

_, _, k_dfc_Standard, k_dfc_adc_Standard = get_kspace(hoseq_Standard; Δt=1)  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_Stitched, k_dfc_adc_Stitched = get_kspace(hoseq_Stitched; Δt=1)  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]


samples_Standard = get_samples(hoseq_Standard; off_val=0) 
samples_Stitched = get_samples(hoseq_Stitched; off_val=0) 
adc_times = KomaMRIBase.get_adc_sampling_times(hoseq_Stitched.SEQ)
t_adc = samples_Stitched.h0.t .* 1e3  # convert to ms  

nSamplePerSegment = Int64(size(t_adc)[1]//4)
nSegment = 4

Segment_range = []
for seg = 1 : nSegment
    r_start = 1+nSamplePerSegment*(seg-1)
    r_stop = seg==nSegment ? nSamplePerSegment*seg : 1+nSamplePerSegment*seg 
    push!(Segment_range, r_start:r_stop)
end


matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 18/2.54
figure_height      = 5.5/2.54
linewidth          = 0.3
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_segment      = ["C0", "C1", "C2", "C3"]

fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(figure_width, figure_height), 
                        facecolor=color_facecoler)
ax1, ax2, ax3 = axs # unpack the 3 axes for zeroth-order, first-order, and second-order plots
vmax = 1.2 * maximum(abs.(samples_Stitched.h1.A*1e3))
vmaxk = 1.1 * maximum([maximum(abs.(k_dfc_adc_Stitched[:,2:3])), maximum(abs.(k_dfc_adc_Standard[:,2:3]))])
for ax in axs
    ax.tick_params(axis="both", length=linewidth*5, width=linewidth, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecoler)
end


for ax in [ax1, ax2]  # remove xticks for first two subplots
    ax.set_xlim(0, t_adc[end]+1)
    ax.set_ylim(-vmax, vmax)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax/4, sigdigits=1)))
    ax.set_xlabel(         "Time (ms)", fontsize=fontsize_label, color=color_label)
    ax.set_ylabel(  "Amplitude (mT/m)", fontsize=fontsize_label, color=color_label)
end
ax3.set_xlabel(        L"k_x (m⁻¹)", fontsize=fontsize_label, color=color_label)
ax3.set_ylabel(        L"k_y (m⁻¹)", fontsize=fontsize_label, color=color_label)

ax1.set_title(L"G_x", fontsize=fontsize_label, color=color_label)
ax2.set_title(L"G_y", fontsize=fontsize_label, color=color_label)
ax3.set_title("k-space trajectory", fontsize=fontsize_label, color=color_label)
ax3.set_xlim(-vmaxk, vmaxk)
ax3.set_ylim(-vmaxk, vmaxk)

for seg = 1:nSegment
    seg_r = Segment_range[seg]
    ax1.plot(t_adc[seg_r], samples_Stitched.h1.A[seg_r]*1e3, color=color_segment[seg], linewidth=linewidth, label="Stitched segment $(seg)")
    ax2.plot(t_adc[seg_r], samples_Stitched.h2.A[seg_r]*1e3, color=color_segment[seg], linewidth=linewidth, label="Stitched segment $(seg)")
end
ax1.plot(t_adc, (samples_Stitched.h1.A-samples_Standard.h1.A)*1e3, color="black", linewidth=linewidth, label="Difference")
ax2.plot(t_adc, (samples_Stitched.h2.A-samples_Standard.h2.A)*1e3, color="black", linewidth=linewidth, label="Difference")



ax3.plot(k_dfc_adc_Stitched[:,2], k_dfc_adc_Stitched[:,3], color="C1", linewidth=linewidth, label="Stitched measurement")
ax3.plot(k_dfc_adc_Standard[:,2], k_dfc_adc_Standard[:,3], color="C0", linewidth=linewidth, label="Single measurement")
ax3.set_aspect(1)

# fig.align_ylabels()
for (ax, ncol) in zip(axs, [1, 1, 1])   # setup legend for each subplot
    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
    loc="upper left", 
    bbox_to_anchor=(0,1.05),
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end
fig.tight_layout(pad=0)
fig.savefig("Figures/Fig3/Fig3.png", dpi=300, bbox_inches="tight")

