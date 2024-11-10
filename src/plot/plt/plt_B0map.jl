"""
    fig = plt_B0map(img::Array{<:Real, 2}; kwargs...)

# Description
    Plots an B0map using matplotlib 

# Arguments
- `img`: (`::Array{<:Real, 2}`) `[nX, nY]`, image array to be plotted

# Keywords
- `title`: (`::Integer`, `=""`) figure's suptitle
- `width`: (`::Real`, `=6`) figure's width
- `height`: (`::Real`, `=5`) figure's height
- `vmax`: (`::Real`, `=100`) maximum value to be used for window width/ window level
- `vmin`: (`::Real`, `=0`) minimum value to be used for window width/ window level
- `cmap`: (`::String`, `="gray"`) colormap to be used for plotting
- `fontsize_title`: (`::Integer`, `=10`) font size of the title
- `fontsize_ticklabel`: (`::Integer`, `=8`) font size of the tick labels
- `pad_labeltick`: (`::Integer`, `=2`) padding between tick labels and the colorbar
- `linewidth`: (`::Real`, `=0.5`) width of the colorbar's outline
- `ticklength`: (`::Real`, `=1.5`) length of the colorbar's ticks
- `color_facecolor`: (`::String`, `="#ffffff"`) background color of the figure

# Returns
- `Figure`: a PyObject representing the figure

# Examples
```julia-repl
julia> fig = plt_B0map(rand(100, 100))

julia> fig.savefig("B0map.png",bbox_inches="tight", pad_inches=0, transparent=true)
```
"""
function plt_B0map(
    b0map::AbstractArray{<:Real, 2};
    title              = ""       ,
    width              = 6        ,
    height             = 5        ,
    vmax               = 100      ,
    vmin               = -100     ,
    cmap               = "jet"    ,
    fontsize_title     = 10       ,
    fontsize_ticklabel = 8        ,
    pad_labeltick      = 2        ,
    linewidth          = 0.5      ,
    ticklength         = 1.5      ,
    color_facecolor    = "#ffffff",
    color_label        = "#cccccc",
    )

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width/2.53999863, height/2.53999863), facecolor=color_facecolor)
    fig.suptitle(title, fontsize=fontsize_title)
    ax.axis("off")
    ai = ax.imshow(b0map, cmap=cmap, vmin=vmin, vmax=vmax)
    cb = fig.colorbar(ai)
    cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
    cb.outline.set_visible(false)
    cb.set_label("B₀ [Hz]", color=color_label, size=fontsize_ticklabel)
    cb.update_ticks()
    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    return fig
end