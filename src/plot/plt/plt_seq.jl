

"""
    p = plt_seq(hoseq::HO_Sequence; kwargs...)

Plots a high order sequence struct.

# Arguments
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct

# Keywords
- `xlabel`: the label for the x-axis.
- `ylabel`: the label for the y-axis.
- `width`: the width of the figure in inches.
- `height`: the height of the figure in inches.
- `fontsize_label`: the font size for the labels.
- `fontsize_legend`: the font size for the legend.
- `fontsize_ticklabel`: the font size for the tick labels.
- `color_facecolor`: the background color of the figure.
- `color_label`: the color of the labels.
- `linewidth`: the width of the curves.
- `ticklength`: the length of the ticks.
- `pad_label`: the padding for the labels.
- `pad_labeltick`: the padding for the tick labels.

# Returns
- `Figure`: a PyObject representing the figure

# Examples
```julia-repl
julia> fig = plt_seq(hoseq)
"""
function plt_seq(
      hoseq ::HO_Sequence             ;
      xlabel             = "Time [ms]",
      ylabel             = "Amplitude",
      width              = 20         ,      
      height             = 8          ,  
      fontsize_label     = 7          ,
      fontsize_legend    = 7          ,
      fontsize_ticklabel = 6          ,
      color_facecolor    = "#1F1F1F"  ,
      color_label        = "#CCCCCC"  ,
      linewidth          = 1.0        ,
      ticklength         = 1.5        ,
      pad_label          = 2          ,
      pad_labeltick      = 2          ,
    )
    sh = SphericalHarmonics()
    # Get the samples of the events in the sequence
    samples = get_samples(hoseq; off_val=Inf, max_rf_samples=Inf)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width/2.53999863, height/2.53999863),facecolor=color_facecolor, squeeze=false)
    ax = ax[1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        # ax.spines[spine].set_linewidth(linewidth)
        ax.spines[spine].set_visible(false)
    end

    # For GRADs
    ax.plot(samples.gx.t*1e3, samples.gx.A*1e3, label="Gx", color="#636EFA", linewidth=linewidth)
    ax.plot(samples.gy.t*1e3, samples.gy.A*1e3, label="Gy", color="#EF553B", linewidth=linewidth)
    ax.plot(samples.gz.t*1e3, samples.gz.A*1e3, label="Gz", color="#00CC96", linewidth=linewidth)


    # For RFs
    O = size(hoseq.SEQ.RF, 1)
    for j in 1:O
        rf_phase = angle.(samples.rf.A[:,j])
        rf_phase[samples.rf.A[:,j] .== Inf] .= Inf
        ax.plot(samples.rf.t*1e3, abs.(samples.rf.A[:,1])*1e6, label="|B1|", color="#AB63FA", linewidth=linewidth)
        ax.plot(samples.rf.t*1e3,                    rf_phase,  label="<B1", color="#FFA15A", linewidth=linewidth, alpha=0.5)
    end

    # For ADCs
    ax.plot(samples.adc.t*1e3, samples.adc.A*1.0, label="ADC", color="#19D3F3", linewidth=linewidth)

    # for dfc measured gradients
    ax.plot(samples.h0.t*1e3, samples.h0.A*1e3, label="h0", color=sh.h0.color, linewidth=linewidth)
    ax.plot(samples.h1.t*1e3, samples.h1.A*1e3, label="h1", color=sh.h1.color, linewidth=linewidth)
    ax.plot(samples.h2.t*1e3, samples.h2.A*1e3, label="h2", color=sh.h2.color, linewidth=linewidth)
    ax.plot(samples.h3.t*1e3, samples.h3.A*1e3, label="h3", color=sh.h3.color, linewidth=linewidth)
    ax.plot(samples.h4.t*1e3, samples.h4.A*1e3, label="h4", color=sh.h4.color, linewidth=linewidth)
    ax.plot(samples.h5.t*1e3, samples.h5.A*1e3, label="h5", color=sh.h5.color, linewidth=linewidth)
    ax.plot(samples.h6.t*1e3, samples.h6.A*1e3, label="h6", color=sh.h6.color, linewidth=linewidth)
    ax.plot(samples.h7.t*1e3, samples.h7.A*1e3, label="h7", color=sh.h7.color, linewidth=linewidth)
    ax.plot(samples.h8.t*1e3, samples.h8.A*1e3, label="h8", color=sh.h8.color, linewidth=linewidth)
    ax.plot(samples.h9.t*1e3, samples.h9.A*1e3, label="h9", color=sh.h9.color, linewidth=linewidth)
    ax.plot(samples.h10.t*1e3, samples.h10.A*1e3, label="h10", color=sh.h10.color, linewidth=linewidth)
    ax.plot(samples.h11.t*1e3, samples.h11.A*1e3, label="h11", color=sh.h11.color, linewidth=linewidth)
    ax.plot(samples.h12.t*1e3, samples.h12.A*1e3, label="h12", color=sh.h12.color, linewidth=linewidth)
    ax.plot(samples.h13.t*1e3, samples.h13.A*1e3, label="h13", color=sh.h13.color, linewidth=linewidth)
    ax.plot(samples.h14.t*1e3, samples.h14.A*1e3, label="h14", color=sh.h14.color, linewidth=linewidth)
    ax.plot(samples.h15.t*1e3, samples.h15.A*1e3, label="h15", color=sh.h15.color, linewidth=linewidth)

    ax.set_xlabel(xlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.set_ylabel(ylabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)

    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, scatteryoffsets=[0.5],
    loc="upper left", bbox_to_anchor=(0.0, 1.00), ncols=15, 
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)
    fig.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
    return fig
end
