using PyPlot
using KomaHighOrder

key = :Standard
hoseq = demo_hoseq(dfc_method=key)   # :Standard or :Stitched
sh = SphericalHarmonics()

_, _, k_dfc, k_dfc_adc = get_kspace(hoseq; Δt=1)  # [nADC, nK]  zeroth-order: [:, 1], first-order: [:, 2:4], and second-order: [:, 5:9]
t_adc = get_adc_sampling_times(hoseq.SEQ) * 1e3 # seconds to milliseconds
t_adc = t_adc .- t_adc[1]

figure_width       = 9
figure_height      = 6
linewidth          = 1
fontsize_legend    = 10
fontsize_label     = 11
fontsize_ticklabel = 8
color_facecolor    = "#ffffff"
color_label        = "#000000"

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(figure_width, figure_height), 
                        sharex=true, facecolor=color_facecolor)
ax1, ax2, ax3 = axs # unpack the 3 axes for zeroth-order, first-order, and second-order plots
vmaxs = 1.2 * [maximum(abs.(k_dfc_adc[:,   1]))*2π, maximum(abs.(k_dfc_adc[:,   2:4])), maximum(abs.(k_dfc_adc[:,   5:9]))]
for (ax, vmax) in zip(axs, vmaxs)
    ax.tick_params(axis="both", color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
    ax.set_xlim(0, t_adc[end]+1)
    ax.set_ylim(-vmax, vmax)
    ax.yaxis.set_major_locator(plt.MultipleLocator(round(vmax/3, sigdigits=1)))
end

for ax in [ax1, ax2]  # remove xticks for first two subplots
    ax.tick_params(axis="x", length=0)
end

ax1.set_ylabel("zeroth-order (rad)", fontsize=fontsize_label, color=color_label)
ax2.set_ylabel( "first-order (m⁻¹)", fontsize=fontsize_label, color=color_label)
ax3.set_ylabel("second-order (m⁻²)", fontsize=fontsize_label, color=color_label)
ax3.set_xlabel(         "time (ms)", fontsize=fontsize_label, color=color_label)

ax1.plot(t_adc, 2π*k_dfc_adc[:, 1], color=sh.h0.color, linewidth=linewidth, label=L"k_{0}")
ax2.plot(t_adc,    k_dfc_adc[:, 2], color=sh.h1.color, linewidth=linewidth, label=L"k_{x}")
ax2.plot(t_adc,    k_dfc_adc[:, 3], color=sh.h2.color, linewidth=linewidth, label=L"k_{y}")
ax2.plot(t_adc,    k_dfc_adc[:, 4], color=sh.h3.color, linewidth=linewidth, label=L"k_{z}")
ax3.plot(t_adc,    k_dfc_adc[:, 5], color=sh.h4.color, linewidth=linewidth, label=L"k_{xy}")
ax3.plot(t_adc,    k_dfc_adc[:, 6], color=sh.h5.color, linewidth=linewidth, label=L"k_{zy}")
ax3.plot(t_adc,    k_dfc_adc[:, 7], color=sh.h6.color, linewidth=linewidth, label=L"k_{2z^2-(x^2+y^2)}")
ax3.plot(t_adc,    k_dfc_adc[:, 8], color=sh.h7.color, linewidth=linewidth, label=L"k_{xz}")
ax3.plot(t_adc,    k_dfc_adc[:, 9], color=sh.h8.color, linewidth=linewidth, label=L"k_{x^2-y^2}")

fig.align_ylabels()
for (ax, ncol) in zip(axs, [1, 3, 5])   # setup legend for each subplot
    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, ncols=ncol, 
    loc="upper left", frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)
end

fig.tight_layout()
fig.savefig("Figures/out/traj_$(String(key)).png", dpi=300)

