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

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width, height), facecolor=color_facecolor)
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