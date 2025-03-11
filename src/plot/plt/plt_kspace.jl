

"""
    p = plt_kspace(hoseq::HO_Sequence; kwargs...)

Plots kspace trajectory of high order sequence struct.

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
julia> fig = plt_kspace(hoseq)
"""
function plt_kspace(
    hoseq ::HO_Sequence             ;
    xlabel             = "kx [1/m]" ,
    ylabel             = "ky [1/m]" ,
    zlabel             = "kz [1/m]" ,
    width              = 20         ,      
    height             = 20         ,  
    fontsize_label     = 7          ,
    fontsize_ticklabel = 6          ,
    color_facecolor    = "#1F1F1F"  ,
    color_label        = "#CCCCCC"  ,
    color_traj         = "#555555"  ,
    linewidth          = 1.0        ,
    ticklength         = 1.5        ,
    pad_label          = 2          ,
    pad_labeltick      = 2          ,
    markersize         = 4          ,
    cmap               = "jet"      ,
    )
    # Get the samples of the events in the sequence
    seq = hoseq.SEQ
	#Calculations of theoretical k-space
	K_nominal, K_nominal_adc, K_dfc, K_dfc_adc = get_kspace(hoseq; Δt=1) #sim_params["Δt"])
    nPoint, nTerm = size(K_nominal_adc)
	K_dfc = K_dfc[:, 2:4]          # H1, H2, H3 => x, y, z
	K_dfc_adc = K_dfc_adc[:, 2:4]

    #Layout
	mink = minimum(K_nominal_adc,dims=1)
	maxk = maximum(K_nominal_adc,dims=1)
	dW = maximum(maxk .- mink, dims=2) * .3
	mink .-= dW
	maxk .+= dW

	t_adc = KomaMRIBase.get_adc_sampling_times(seq)

	fig = plt.figure(figsize=(width/2.53999863, height/2.53999863),facecolor=color_facecolor)
	ax = fig.add_subplot(1, 1, 1, projection="3d")
    ax.set_box_aspect((1,1,1)) 
    ax.grid(false) 
    ax.xaxis.set_pane_color(color_facecolor)
    ax.yaxis.set_pane_color(color_facecolor)
    ax.zaxis.set_pane_color(color_facecolor)
    ax.xaxis.line.set_color(color_label) 
    ax.yaxis.line.set_color(color_label) 
    ax.zaxis.line.set_color(color_label)
    ax.set_xlim(mink[1], maxk[1])
    ax.set_ylim(mink[2], maxk[2])
    ax.set_zlim(mink[3], maxk[3])
    ax.xaxis.set_tick_params("both", colors=color_label, length=ticklength, width=linewidth, pad=pad_labeltick, 
    color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    ax.yaxis.set_tick_params("both", colors=color_label, length=ticklength, width=linewidth, pad=pad_labeltick, 
    color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    ax.zaxis.set_tick_params("both", colors=color_label, length=ticklength, width=linewidth, pad=pad_labeltick, 
    color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)

    ax.set_facecolor(color_facecolor)
    # ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
    #     color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    # for spine in ax.spines  # "left", "right", "bottom", "top"
    #     # ax.spines[spine].set_linewidth(linewidth)
    #     ax.spines[spine].set_visible(false)
    # end
    ax.plot(K_nominal[:,1], K_nominal[:,2], K_nominal[:,3], color=color_traj, linewidth=linewidth)
    ax.scatter(K_nominal_adc[:,1], K_nominal_adc[:,2], K_nominal_adc[:,3], marker=".", s=markersize, c=collect(1:nPoint), cmap=cmap)

    ax.set_xlabel(xlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.set_ylabel(ylabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.set_zlabel(zlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    fig.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
    return fig
end

