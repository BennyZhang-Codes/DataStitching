
"""
    fig = plt_image(img::Array{<:Real, 2}; title="", width=5, height=5, vmaxp=100, vminp=0, cmap="gray", fontsize_title=10, color_facecolor="#ffffff")
    Plots an image array using matplotlib 

# Arguments
- `img`: (`::Array{<:Real, 2}`) `[nX, nY]`, image array to be plotted

# Keywords
- `title`: (`::Integer`, `=""`) figure's suptitle
- `width`: (`::Real`, `=5`) figure's width
- `height`: (`::Real`, `=5`) figure's height
- `vmaxp`: (`::Real`, `=100`) percentile of maximum value to be used for window width/ window level
- `vminp`: (`::Real`, `=0`) percentile of minimum value to be used for window width/ window level
- `cmap`: (`::String`, `="gray"`) colormap to be used for plotting
- `fontsize_title`: (`::Integer`, `=10`) font size of the title
- `color_facecolor`: (`::String`, `="#ffffff"`) background color of the figure

# Returns
- `Figure`: a PyObject representing the figure

# Examples
```julia-repl
julia> fig = plt_image(rand(100, 100))

julia> fig.savefig("123.png",bbox_inches="tight", pad_inches=0, transparent=true)
```
"""
function plt_image(
    img::AbstractArray{<:Real, 2};
    title              = ""       ,
    width              = 5        ,
    height             = 5        ,
    vmaxp              = 100      ,
    vminp              = 0        ,
    cmap               = "gray"   ,
    fontsize_title     = 10       ,
    color_facecolor    = "#ffffff",
    )
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width, height), facecolor=color_facecolor)
	fig.suptitle(title, fontsize=fontsize_title)
    ax.imshow(img, cmap=cmap, vmin=quantile(reshape(img, :), vminp/100), vmax=quantile(reshape(img, :), vmaxp/100))
    ax.axis("off")
    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    return fig
end

"""
    fig = plt_images(imgs::Array{<:Real, 3}; title="", width=5, height=5, vmaxp=100, vminp=0, cmap="gray", fontsize_title=10, color_facecolor="#ffffff")
    Plots a sequence of image array using matplotlib 

# Arguments
- `imgs`: (`::Array{<:Real, 3}`) `[nFrame, nX, nY]`, image array to be plotted

# Keywords
- `title`: (`::String`, `=""`) figure's suptitle
- `width`: (`::Real`, `=5`) figure's width
- `height`: (`::Real`, `=5`) figure's height
- `vmaxp`: (`::Real`, `=100`) percentile of maximum value to be used for window width/ window level
- `vminp`: (`::Real`, `=0`) percentile of minimum value to be used for window width/ window level
- `cmap`: (`::String`, `="gray"`) colormap to be used for plotting
- `fontsize_title`: (`::Integer`, `=10`) font size of the title
- `color_facecolor`: (`::String`, `="#ffffff"`) background color of the figure

# Returns
- `Figure`: a PyObject representing the figure

# Examples
```julia-repl
julia> imgs = rand(10, 100, 100)

julia> fig = plt_images(imgs)

julia> fig.savefig("123.png",bbox_inches="tight", pad_inches=0, transparent=true)
"""
function plt_images(
    imgs::AbstractArray{<:Real, 3};
    nRow               = nothing  ,
    nCol               = nothing  ,
    title              = ""       ,
    width              = 5        ,
    height             = 5        ,
    vmaxp              = 100      ,
    vminp              = 0        ,
    cmap               = "gray"   ,
    fontsize_title     = 10       ,
    color_facecolor    = "#ffffff",
    )
    nFrame, nX, nY = size(imgs)
    if nRow === nothing || nCol === nothing
        nRow, nCol = get_factors(nFrame)
    end

    vmin, vmax = quantile(reshape(imgs, :), vminp/100), quantile(reshape(imgs, :), vmaxp/100)

    fig, axs = plt.subplots(nrows=nRow, ncols=nCol, figsize=(width, height), facecolor=color_facecolor)
	axs = reshape(axs, nRow, nCol) # to keep axs as a 2D array
    
    fig.suptitle(title, fontsize=fontsize_title)
    for frame = 1:nFrame
        row, col = Int(ceil(frame/nCol)), (frame-1)%nCol + 1
        ax = axs[row, col]
        ax.imshow(imgs[frame, :, :], cmap=cmap, vmin=vmin, vmax=vmax)
        ax.axis("off")
    end
    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    return fig
end



