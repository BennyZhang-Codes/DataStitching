

"""
    p = plt_seqd(hoseq::HO_Sequence; kwargs...)

Plots a high order discretized sequence struct.

# Arguments
- `hoseqd`: (`::HHO_DiscreteSequence `) HO_DiscreteSequence  struct
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct

# Keywords
- `sampling_params`: The sampling parameters for the sequence.
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
- `marker_size`: the size of the markers.

# Returns
- `Figure`: a PyObject representing the figure

# Examples
```julia-repl
julia> fig = plt_seqd(hoseq)
"""
function plt_seqd(
    hoseqd :: HO_DiscreteSequence                             ;
    xlabel             = "Time [ms]"                          ,
    ylabel             = "Amplitude"                          ,
    width              = 20                                   ,      
    height             = 8                                    ,  
    fontsize_label     = 7                                    ,
    fontsize_legend    = 7                                    ,
    fontsize_ticklabel = 6                                    ,
    color_facecolor    = "#1F1F1F"                            ,
    color_label        = "#CCCCCC"                            ,
    linewidth          = 1.0                                  ,
    ticklength         = 1.5                                  ,
    pad_label          = 2                                    ,
    pad_labeltick      = 2                                    ,
    marker_size        = 4                                    ,
    )
    sh = SphericalHarmonics()

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width/2.53999863, height/2.53999863),facecolor=color_facecolor, squeeze=false)
    ax = ax[1]
    ax.set_facecolor(color_facecolor)
    ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        # ax.spines[spine].set_linewidth(linewidth)
        ax.spines[spine].set_visible(false)
    end

    ax.plot(hoseqd.seqd.t*1e3, hoseqd.seqd.Gx*1e3, label="Gx", marker="o", color="#636efa", markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.seqd.Gy*1e3, label="Gy", marker="o", color="#EF553B", markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.seqd.Gz*1e3, label="Gz", marker="o", color="#00cc96", markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, abs.(hoseqd.seqd.B1*1e6), label="|B1|", marker="o", color="#AB63FA", markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t[hoseqd.seqd.ADC]*1e3, zeros(sum(hoseqd.seqd.ADC)), label="ADC", marker="x", color="#FF6692", markersize=marker_size, linewidth=linewidth)

    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h0*1e3, label="h0", marker="o", color=sh.h0.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h1*1e3, label="h1", marker="o", color=sh.h1.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h2*1e3, label="h2", marker="o", color=sh.h2.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h3*1e3, label="h3", marker="o", color=sh.h3.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h4*1e3, label="h4", marker="o", color=sh.h4.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h5*1e3, label="h5", marker="o", color=sh.h5.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h6*1e3, label="h6", marker="o", color=sh.h6.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h7*1e3, label="h7", marker="o", color=sh.h7.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h8*1e3, label="h8", marker="o", color=sh.h8.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h9*1e3, label="h9", marker="o", color=sh.h9.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h10*1e3, label="h10", marker="o", color=sh.h10.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h11*1e3, label="h11", marker="o", color=sh.h11.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h12*1e3, label="h12", marker="o", color=sh.h12.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h13*1e3, label="h13", marker="o", color=sh.h13.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h14*1e3, label="h14", marker="o", color=sh.h14.color, markersize=marker_size, linewidth=linewidth)
    ax.plot(hoseqd.seqd.t*1e3, hoseqd.h15*1e3, label="h15", marker="o", color=sh.h15.color, markersize=marker_size, linewidth=linewidth)

    ax.set_xlabel(xlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.set_ylabel(ylabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)

    ax.legend(fontsize=fontsize_legend, labelcolor=color_label, scatteryoffsets=[0.5],
    loc="upper left", bbox_to_anchor=(0.0, 1.00), ncols=15, 
    frameon=false, handlelength=1, handletextpad=0.5, columnspacing=1)
    fig.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
    return fig
end



function plt_seqd(
    hoseq ::HO_Sequence                                       ;
    sampling_params    = KomaMRIBase.default_sampling_params(),
    xlabel             = "Time [ms]"                          ,
    ylabel             = "Amplitude"                          ,
    width              = 20                                   ,      
    height             = 8                                    ,  
    fontsize_label     = 7                                    ,
    fontsize_legend    = 7                                    ,
    fontsize_ticklabel = 6                                    ,
    color_facecolor    = "#1F1F1F"                            ,
    color_label        = "#CCCCCC"                            ,
    linewidth          = 1.0                                  ,
    ticklength         = 1.5                                  ,
    pad_label          = 2                                    ,
    pad_labeltick      = 2                                    ,
    marker_size        = 4                                    ,
    )
    hoseqd = discretize(hoseq; sampling_params)
    return plt_seqd(hoseqd; xlabel, ylabel, width, height, fontsize_label, fontsize_legend, fontsize_ticklabel, 
    color_facecolor, color_label, linewidth, ticklength, pad_label, pad_labeltick, marker_size)
end