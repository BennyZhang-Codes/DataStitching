using PyPlot
using KomaHighOrder
using LinearAlgebra

sh = SphericalHarmonics()
nSegment_fullysampled = 4
nSegment_undersampled = 36

_, ntStitched_fullysampled, _ = load_dfc(;dfc_method=:Stitched, seqname="demo", r=1)
_, ntStitched_undersampled, _ = load_dfc(;dfc_method=:Stitched, seqname="demo", r=30)

# seq
hoseq_fullysampled_Standard = demo_hoseq(dfc_method=:Standard)[8]   # :Standard or :Stitched
hoseq_fullysampled_Stitched = demo_hoseq(dfc_method=:Stitched)[8]   # :Standard or :Stitched
hoseq_undersampled_Standard = demo_hoseq(dfc_method=:Standard, r=30)[8]   # :Standard or :Stitched
hoseq_undersampled_Stitched = demo_hoseq(dfc_method=:Stitched, r=30)[8]   # :Standard or :Stitched

# gradient
samples_fullysampled_Standard = get_samples(hoseq_fullysampled_Standard; off_val=0); 
samples_fullysampled_Stitched = get_samples(hoseq_fullysampled_Stitched; off_val=0); 
samples_undersampled_Standard = get_samples(hoseq_undersampled_Standard; off_val=0); 
samples_undersampled_Stitched = get_samples(hoseq_undersampled_Stitched; off_val=0); 

adc_times_fullysampled = KomaMRIBase.get_adc_sampling_times(hoseq_fullysampled_Stitched.SEQ);
adc_times_undersampled = KomaMRIBase.get_adc_sampling_times(hoseq_undersampled_Stitched.SEQ);
t_adc_fullysampled = samples_fullysampled_Stitched.h0.t .* 1e3;  # convert to ms  
t_adc_undersampled = samples_undersampled_Stitched.h0.t .* 1e3;  # convert to ms  

# trajectory
_, _, k_dfc_fullysampled_Standard, k_dfc_adc_fullysampled_Standard = get_kspace(hoseq_fullysampled_Standard; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_fullysampled_Stitched, k_dfc_adc_fullysampled_Stitched = get_kspace(hoseq_fullysampled_Stitched; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_undersampled_Standard, k_dfc_adc_undersampled_Standard = get_kspace(hoseq_undersampled_Standard; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_undersampled_Stitched, k_dfc_adc_undersampled_Stitched = get_kspace(hoseq_undersampled_Stitched; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]


# nSamplePerSegment_fullysampled = Int64(size(t_adc_fullysampled)[1]//4)
# Segment_range_fullysampled = []
# for seg = 1 : nSegment_fullysampled
#     r_start = 1+nSamplePerSegment_fullysampled*(seg-1)
#     r_stop = seg==nSegment_fullysampled ? nSamplePerSegment_fullysampled*seg : 1+nSamplePerSegment_fullysampled*seg 
#     push!(Segment_range_fullysampled, r_start:r_stop)
# end
# nSamplePerSegment_undersampled = Int64(size(t_adc_undersampled)[1]//36)
# Segment_range_undersampled = []
# for seg = 1 : nSegment_undersampled
#     r_start = 1+nSamplePerSegment_undersampled*(seg-1)
#     r_stop = seg==nSegment_undersampled ? nSamplePerSegment_undersampled*seg : 1+nSamplePerSegment_undersampled*seg 
#     push!(Segment_range_undersampled, r_start:r_stop)
# end

a = Int64.(vec(ntStitched_fullysampled))
a[1] += 3
e = cumsum(a)
s = e .- a
Segment_range_fullysampled = diag([st+1:en+1 for st in s, en in e])


a = Int64.(vec(ntStitched_undersampled))
a[1] += 3
e = cumsum(a)
s = e .- a
Segment_range_undersampled = diag([st+1:en+1 for st in s, en in e])

matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 18/2.54
figure_height      = 11/2.54
linewidth          = 0.4
fontsize_legend    = 5
fontsize_label     = 6
fontsize_ticklabel = 4
color_facecoler    = "#ffffff"
color_label        = "#000000"
color_segment      = ["C0", "C1", "C2", "C3"]

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(figure_width, figure_height), 
                        facecolor=color_facecoler)
ax1, ax2, ax3, ax4, ax5, ax6 = axs # unpack the 3 axes for zeroth-order, first-order, and second-order plots
vmax_g_fullysampled = 1.2 * maximum(abs.(samples_fullysampled_Stitched.h1.A*1e3))
vmax_k_fullysampled = 1.1 * maximum([maximum(abs.(k_dfc_adc_fullysampled_Stitched[:,2:3])), maximum(abs.(k_dfc_adc_fullysampled_Standard[:,2:3]))])
vmax_g_undersampled = 1.2 * maximum(abs.(samples_undersampled_Stitched.h1.A*1e3))
vmax_k_undersampled = 1.1 * maximum([maximum(abs.(k_dfc_adc_undersampled_Stitched[:,2:3])), maximum(abs.(k_dfc_adc_undersampled_Standard[:,2:3]))])

titles = [L"G_x", L"G_x", L"G_y", L"G_y", "k-space trajectory", "k-space trajectory"]
orders = ["a", "d", "b", "e", "c", "f"]
for (ax, title, order) in zip(axs, titles, orders)
    ax.set_title(title, fontsize=fontsize_label, color=color_label)
    ax.text(0, 1.1, "($(order))", fontsize=fontsize_label, color=color_label, transform=ax.transAxes, ha="left", va="top")
end


for ax in axs
    ax.tick_params(axis="both", length=linewidth*5, width=linewidth, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecoler)
end
for ax in axs[1, 1:2] 
    ax.set_xlim(-1, t_adc_fullysampled[end]+1)
    ax.set_ylim(-vmax_g_fullysampled, vmax_g_fullysampled)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax_g_fullysampled/4, sigdigits=1)))
end
for ax in axs[2, 1:2] 
    ax.set_xlim(-1, t_adc_undersampled[end]+1)
    ax.set_ylim(-vmax_g_undersampled, vmax_g_undersampled)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax_g_undersampled/4, sigdigits=1)))
end
for ax in axs[:, 1:2]  # remove xticks for first two subplots
    ax.set_xlabel(         "Time (ms)", fontsize=fontsize_label, color=color_label)
    ax.set_ylabel(  "Amplitude (mT/m)", fontsize=fontsize_label, color=color_label)
end
for ax in axs[:, 3]  # remove yticks for first two subplots
    ax.set_xlabel(L"k_x (m⁻¹)", fontsize=fontsize_label, color=color_label)
    ax.set_ylabel(L"k_y (m⁻¹)", fontsize=fontsize_label, color=color_label)
end

for seg = 1:nSegment_fullysampled
    seg_r = Segment_range_fullysampled[seg]
    axs[1,1].plot(t_adc_fullysampled[seg_r], samples_fullysampled_Stitched.h1.A[seg_r]*1e3, color=color_segment[seg], linewidth=linewidth, label="Stitched segment $(seg)")
    axs[1,2].plot(t_adc_fullysampled[seg_r], samples_fullysampled_Stitched.h2.A[seg_r]*1e3, color=color_segment[seg], linewidth=linewidth, label="Stitched segment $(seg)")
end
axs[1,1].plot(t_adc_fullysampled, (samples_fullysampled_Stitched.h1.A-samples_fullysampled_Standard.h1.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[1,2].plot(t_adc_fullysampled, (samples_fullysampled_Stitched.h2.A-samples_fullysampled_Standard.h2.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[1,3].plot(k_dfc_adc_fullysampled_Stitched[:,2], k_dfc_adc_fullysampled_Stitched[:,3], color="C1", linewidth=linewidth, label="Stitched measurement")
axs[1,3].plot(k_dfc_adc_fullysampled_Standard[:,2], k_dfc_adc_fullysampled_Standard[:,3], color="C0", linewidth=linewidth, label="Single measurement")
axs[1,3].set_aspect(1)

for seg = 1:nSegment_undersampled
    seg_r = Segment_range_undersampled[seg]
    axs[2,1].plot(t_adc_undersampled[seg_r], samples_undersampled_Stitched.h1.A[seg_r]*1e3, linewidth=linewidth)
    axs[2,2].plot(t_adc_undersampled[seg_r], samples_undersampled_Stitched.h2.A[seg_r]*1e3, linewidth=linewidth)
end
axs[2,1].plot(t_adc_undersampled, (samples_undersampled_Stitched.h1.A-samples_undersampled_Standard.h1.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[2,2].plot(t_adc_undersampled, (samples_undersampled_Stitched.h2.A-samples_undersampled_Standard.h2.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[2,3].plot(k_dfc_adc_undersampled_Stitched[:,2], k_dfc_adc_undersampled_Stitched[:,3], color="C1", linewidth=linewidth, label="Stitched measurement")
axs[2,3].plot(k_dfc_adc_undersampled_Standard[:,2], k_dfc_adc_undersampled_Standard[:,3], color="C0", linewidth=linewidth, label="Single measurement")
axs[2,3].set_aspect(1)



# fig.align_ylabels()
for (ax, ncol) in zip(axs, ones(6)*1)   # setup legend for each subplot
    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
    loc="upper left", 
    bbox_to_anchor=(0,1.05),
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end
fig.align_ylabels()
fig.tight_layout(pad=0, h_pad=0.5, w_pad=0.5)
fig.savefig("Figures/Fig3/Fig3.png", dpi=300, bbox_inches="tight")

