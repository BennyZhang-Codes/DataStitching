

"""
    p = plt_magnetization(hoseq::HO_Phantom; kwargs...)

Plots a high order phantom object.

# Arguments
- `ph`: (`::HO_Phantom`) HO_Phantom struct
- `mag`: (`::Mag`) Magnetization struct
- `key`: (`::Symbol`) key for the phantom to be plotted

# Keywords
- `t0`: (`::Real`) time offset for the phantom
- `view_2d`: (`::Bool`) whether to plot the phantom in 2D
- `xlabel`: the label for the x-axis.
- `ylabel`: the label for the y-axis.
- `zlabel`: the label for the z-axis.
- `width`: the width of the figure in inches.
- `height`: the height of the figure in inches.
- `fontsize_label`: the font size for the labels.
- `fontsize_ticklabel`: the font size for the tick labels.
- `color_facecolor`: the background color of the figure.
- `color_label`: the color of the labels.
- `linewidth`: the width of the curves.
- `ticklength`: the length of the ticks.
- `pad_label`: the padding for the labels.
- `pad_labeltick`: the padding for the tick labels.
- `markersize`: the size of the markers.

# Returns
- `Figure`: a PyObject representing the figure

# Examples
```julia-repl
julia> fig = plt_magnetization(phantom, mag, :xy_mag)
"""
function plt_magnetization(
    ph                 :: HO_Phantom,
    mag                :: Mag       ,
    key                :: Symbol    ;
    zmin               = nothing    ,
	zmax               = nothing    ,
    t0                 = 0          ,
    view_2d            = false      ,
    xlabel             = "x [cm]"   ,
    ylabel             = "y [cm]"   ,
    zlabel             = "z [cm]"   ,
    width              = 20         ,      
    height             = 20         ,  
    fontsize_label     = 7          ,
    fontsize_ticklabel = 6          ,
    color_facecolor    = "#1F1F1F"  ,
    color_label        = "#CCCCCC"  ,
    linewidth          = 1.0        ,
    ticklength         = 1.5        ,
    pad_label          = 2          ,
    pad_labeltick      = 2          ,
    markersize         = 4          ,
    )
    @assert key in [:xy_mag, :xy_pha, :xy_real, :xy_imag, :z] "key must be :xy_mag, :xy_pha, :xy_real, xy_imag, or :z"
    @assert size(mag.xy, 1) == size(mag.z, 1) == size(ph.x, 1) "Mag and Phantom must have the same number of spins."

	if key == :xy_mag
		spinstate = abs.(mag.xy)
    elseif key == :xy_pha
        spinstate = angle.(mag.xy)
    elseif key == :xy_real
        spinstate = real.(mag.xy)
    elseif key == :xy_imag
        spinstate = imag.(mag.xy)
    else
        spinstate = mag.z
    end

	cmin_key = isnothing(zmin) ? minimum(spinstate) : zmin
	cmax_key = isnothing(zmax) ? maximum(spinstate) : zmax
	unit   = ""
	cmap   = "gray"
    factor = 1.0

	x0 = -maximum(abs.([ph.x ph.y ph.z]))*1e2
    xf =  maximum(abs.([ph.x ph.y ph.z]))*1e2


	fig = plt.figure(figsize=(width/2.53999863, height/2.53999863),facecolor=color_facecolor)

    if view_2d
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect("equal", "box")

        ax.set_facecolor(color_facecolor)
        ax.tick_params(axis="both", length=ticklength, width=linewidth, pad=pad_labeltick, 
            color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
        for spine in ax.spines  # "left", "right", "bottom", "top"
            # ax.spines[spine].set_linewidth(linewidth)
            ax.spines[spine].set_visible(false)
        end
        ax.set_xlim(x0,xf)
        ax.set_ylim(x0,xf)
        
        sc = ax.scatter(
            (ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
            (ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
            marker=".", s=markersize, c=spinstate*factor,  
            cmap=cmap, vmin=cmin_key, vmax=cmax_key
        )
        cb = fig.colorbar(sc, ax=ax, shrink=0.6)
        cb.ax.set_title("[$(String(key))]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
        cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
        cb.outline.set_visible(false)
        cb.update_ticks()

    else
        ax = fig.add_subplot(1, 1, 1, projection="3d")
        ax.set_box_aspect((1,1,1)) 
        ax.grid(false) 
        ax.xaxis.set_pane_color(color_facecolor)
        ax.yaxis.set_pane_color(color_facecolor)
        ax.zaxis.set_pane_color(color_facecolor)
        ax.xaxis.line.set_color(color_label) 
        ax.yaxis.line.set_color(color_label) 
        ax.zaxis.line.set_color(color_label)
        ax.set_xlim(x0,xf)
        ax.set_ylim(x0,xf)
        ax.set_zlim(x0,xf)

        ax.xaxis.set_tick_params("both", colors=color_label, length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
        ax.yaxis.set_tick_params("both", colors=color_label, length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
        ax.zaxis.set_tick_params("both", colors=color_label, length=ticklength, width=linewidth, pad=pad_labeltick, 
        color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
    
        ax.set_facecolor(color_facecolor)

        sc = ax.scatter(
            (ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
            (ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
            (ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2, 
            marker=".", s=markersize, c=spinstate*factor,  
            cmap=cmap, vmin=cmin_key, vmax=cmax_key
        )
        cb = fig.colorbar(sc, ax=ax, shrink=0.6)
        cb.ax.set_title("[$(String(key))]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
        cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
        cb.outline.set_visible(false)
        cb.update_ticks()

        ax.set_zlabel(zlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    end

    ax.set_xlabel(xlabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    ax.set_ylabel(ylabel, fontsize=fontsize_label, color=color_label, labelpad=pad_label)
    
    fig.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
    return fig
end
