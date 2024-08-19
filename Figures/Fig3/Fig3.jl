using PyPlot
using KomaHighOrder
using LinearAlgebra

sh = SphericalHarmonics()
nSegment_fully = 4
nSegment_under = 36

_, ntStitched_fully, _ = load_dfc(;dfc_method=:Stitched, seqname="demo", r=1)
_, ntStitched_under, _ = load_dfc(;dfc_method=:Stitched, seqname="demo", r=30)

# seq
hoseq_fully_Standard = demo_hoseq(dfc_method=:Standard)[8]   # :Standard or :Stitched
hoseq_fully_Stitched = demo_hoseq(dfc_method=:Stitched)[8]   # :Standard or :Stitched
hoseq_under_Standard = demo_hoseq(dfc_method=:Standard, r=30)[8]   # :Standard or :Stitched
hoseq_under_Stitched = demo_hoseq(dfc_method=:Stitched, r=30)[8]   # :Standard or :Stitched

# gradient
samples_fully_Standard = get_samples(hoseq_fully_Standard; off_val=0); 
samples_fully_Stitched = get_samples(hoseq_fully_Stitched; off_val=0); 
samples_under_Standard = get_samples(hoseq_under_Standard; off_val=0); 
samples_under_Stitched = get_samples(hoseq_under_Stitched; off_val=0); 

adc_times_fully = KomaMRIBase.get_adc_sampling_times(hoseq_fully_Stitched.SEQ);
adc_times_under = KomaMRIBase.get_adc_sampling_times(hoseq_under_Stitched.SEQ);
t_adc_fully = samples_fully_Stitched.h0.t .* 1e3;  # convert to ms  
t_adc_under = samples_under_Stitched.h0.t .* 1e3;  # convert to ms  

# trajectory
_, _, k_dfc_fully_Standard, k_dfc_adc_fully_Standard = get_kspace(hoseq_fully_Standard; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_fully_Stitched, k_dfc_adc_fully_Stitched = get_kspace(hoseq_fully_Stitched; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_under_Standard, k_dfc_adc_under_Standard = get_kspace(hoseq_under_Standard; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
_, _, k_dfc_under_Stitched, k_dfc_adc_under_Stitched = get_kspace(hoseq_under_Stitched; Δt=1);  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]


# nSamplePerSegment_fully = Int64(size(t_adc_fully)[1]//4)
# Segment_range_fully = []
# for seg = 1 : nSegment_fully
#     r_start = 1+nSamplePerSegment_fully*(seg-1)
#     r_stop = seg==nSegment_fully ? nSamplePerSegment_fully*seg : 1+nSamplePerSegment_fully*seg 
#     push!(Segment_range_fully, r_start:r_stop)
# end

a = Int64.(vec(ntStitched_fully))
a[1] += 3
e = cumsum(a)
s = e .- a
Segment_range_fully = diag([st+1:en+1 for st in s, en in e])


a = Int64.(vec(ntStitched_under))
a[1] += 3
e = cumsum(a)
s = e .- a
Segment_range_under = diag([st+1:en+1 for st in s, en in e])

matplotlib.rc("mathtext", default="regular")
matplotlib.rc("figure", dpi=200)
matplotlib.rc("font", family="Times New Roman")
matplotlib.rcParams["mathtext.default"]
figure_width       = 18/2.54
figure_height      = 11/2.54
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
color_segment      = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",  "#e377c2",  "#bcbd22","#17becf"]#,"#8c564b"  ,"#7f7f7f"
fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(figure_width, figure_height), 
                        facecolor=color_facecoler)
ax1, ax2, ax3, ax4, ax5, ax6 = axs # unpack the 3 axes for zeroth-order, first-order, and second-order plots
vmax_g_fully = 1.1 * maximum(abs.(samples_fully_Stitched.h1.A*1e3))
vmax_k_fully = 1.1 * maximum([maximum(abs.(k_dfc_adc_fully_Stitched[:,2:3])), maximum(abs.(k_dfc_adc_fully_Standard[:,2:3]))])
vmax_g_under = 1.1 * maximum(abs.(samples_under_Stitched.h1.A*1e3))
vmax_k_under = 1.1 * maximum([maximum(abs.(k_dfc_adc_under_Stitched[:,2:3])), maximum(abs.(k_dfc_adc_under_Standard[:,2:3]))])


for ax in axs
    ax.tick_params(axis="both", length=ticklength, width=linewidth, 
    color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecoler)
end

titles = [L"G_x", L"G_y", "k-space trajectory"]
xlabels = ["Time (ms)", "Time (ms)", L"k_x (m⁻¹)"]
ylabels = ["Amplitude (mT/m)", "Amplitude (mT/m)", L"k_y (m⁻¹)"]
for row = 1:2
    for col = 1:3
        ax = axs[row,col]
        title = titles[col]
        xlabel = xlabels[col]
        ylabel = ylabels[col]
        ax.set_title(title, fontsize=fontsize_label, color=color_label)
        ax.set_xlabel(xlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
        ax.set_ylabel(ylabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    end
end

for ax in axs[1, 1:2] 
    ax.set_xlim(t_adc_fully[1], t_adc_fully[end])
    ax.set_ylim(-vmax_g_fully, vmax_g_fully)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax_g_fully/4, sigdigits=1)))
end
for ax in axs[2, 1:2] 
    ax.set_xlim(t_adc_under[1], t_adc_under[end])
    ax.set_ylim(-vmax_g_under, vmax_g_under)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax_g_under/4, sigdigits=1)))
end

# fully
for seg = 1:nSegment_fully
    seg_r = Segment_range_fully[seg]
    axs[1,1].plot(t_adc_fully[seg_r], samples_fully_Stitched.h1.A[seg_r]*1e3, color=color_segment[seg], linewidth=linewidth, label="Stitched segment $(seg)")
    axs[1,2].plot(t_adc_fully[seg_r], samples_fully_Stitched.h2.A[seg_r]*1e3, color=color_segment[seg], linewidth=linewidth, label="Stitched segment $(seg)")
end
axs[1,1].plot(t_adc_fully, (samples_fully_Stitched.h1.A-samples_fully_Standard.h1.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[1,2].plot(t_adc_fully, (samples_fully_Stitched.h2.A-samples_fully_Standard.h2.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[1,3].plot(k_dfc_adc_fully_Stitched[:,2], k_dfc_adc_fully_Stitched[:,3], color="C1", linewidth=linewidth, label="Stitched measurement")
axs[1,3].plot(k_dfc_adc_fully_Standard[:,2], k_dfc_adc_fully_Standard[:,3], color="C0", linewidth=linewidth, label="Single measurement")
axs[1,3].set_aspect(1)

# under
for seg = 1:nSegment_under
    nColor = length(color_segment)
    color_idx = seg%nColor != 0 ? seg%nColor : nColor
    
    seg_r = Segment_range_under[seg]
    axs[2,1].plot(t_adc_under[seg_r], samples_under_Stitched.h1.A[seg_r]*1e3, color=color_segment[color_idx], linewidth=linewidth, label="Stitched measurement")
    axs[2,2].plot(t_adc_under[seg_r], samples_under_Stitched.h2.A[seg_r]*1e3, color=color_segment[color_idx], linewidth=linewidth, label="Stitched measurement")
end 
lDiff_gx, = axs[2,1].plot(t_adc_under, (samples_under_Stitched.h1.A-samples_under_Standard.h1.A)*1e3, color="black", linewidth=linewidth, label="Difference")
lDiff_gy, = axs[2,2].plot(t_adc_under, (samples_under_Stitched.h2.A-samples_under_Standard.h2.A)*1e3, color="black", linewidth=linewidth, label="Difference")
axs[2,3].plot(k_dfc_adc_under_Stitched[:,2], k_dfc_adc_under_Stitched[:,3], color="C1", linewidth=linewidth, label="Stitched measurement")
axs[2,3].plot(k_dfc_adc_under_Standard[:,2], k_dfc_adc_under_Standard[:,3], color="C0", linewidth=linewidth, label="Single measurement")
axs[2,3].set_aspect(1)


# legend
for (ax, ncol) in zip(axs[[1,3, 5, 6]], [2, 2, 2, 2])   # setup legend for each subplot
    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
    loc="upper left", 
    bbox_to_anchor=(0,1.05),
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end
ls = []
for color in color_segment
    push!(ls, lines.Line2D([], [], color=color, linewidth=linewidth))
end
for (ax, lDiff) in zip(axs[2,1:2], [lDiff_gx, lDiff_gy])
    ax.legend([Tuple(ls), lDiff], ["Stitched measurement", "Difference"] , handler_map=Dict(Tuple(ls) => legend_handler.HandlerTuple(ndivide=length(color_segment), pad=0)),
        fontsize=fontsize_legend, labelcolor=color_label, ncols=2, loc="upper left", bbox_to_anchor=(0,1.05),
        frameon=false, handlelength=0.4*length(color_segment), handletextpad=0.5, columnspacing=1,labelspacing=0.2)
end



fig.align_ylabels()
orders = ["a", "b", "c", "d", "e", "f"]
for row = 1:2
    for col = 1:3
        order = orders[3*(row-1)+col]
        fig.text(0.01+(col-1)/3, 1-0.5*(row-1), "($(order))", ha="left", va="baseline", fontsize=fontsize_subfigure)
    end
end
fig.tight_layout(pad=0, h_pad=0.5, w_pad=0.5)
fig.savefig("Figures/Fig3/Fig3.png", dpi=300, bbox_inches="tight")