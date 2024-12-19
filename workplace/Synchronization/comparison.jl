using PyPlot
using PyCall
animation = pyimport("matplotlib.animation")


img1 = randn(100, 100);
img2 = randn(100, 100);
imgs = [img1, img2]

fig, ax = plt.subplots()
im1 = ax.imshow(imgs[1])

function update(frame)
    ax.set_title("Frame $frame")
    im1.set_data(imgs[frame%2+1])
    return [im1]
end 
ani = animation.FuncAnimation(fig, update, frames=20, interval=1000, blit=true)

ani.save("rotation_animation.gif", writer="imagemagick", fps=30, dpi=300)
# plt.show()

function plt_imagec(
    img1::AbstractArray{<:Real, 2},
    img2::AbstractArray{<:Real, 2};
    nFrames            = 10       ,
    interval           = 50       ,
    title              = ""       ,
    width              = 5        ,
    height             = 5        ,
    vmaxp              = 100      ,
    vminp              = 0        ,
    vmax               = nothing  ,
    vmin               = nothing  ,
    cmap               = "gray"   ,
    fontsize_title     = 10       ,
    color_facecolor    = "#ffffff",
    color_label        = "#000000",
    )
    @assert size(img1) == size(img2) "The size of two images should be the same."
    imgs = [img1, img2]
    nX, nY = size(img1)
    if vmax === nothing || vmin === nothing
        vmax=quantile(reshape([img1; img2], :), vmaxp/100)
        vmin=quantile(reshape([img1; img2], :), vminp/100)
    end

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(width/2.53999863, (nY/nX)*height/2.53999863), facecolor=color_facecolor)
    ax.set_title(title, fontsize=fontsize_title, color=color_label)
    im1 = ax.imshow(imgs[1], cmap=cmap, vmin=vmin, vmax=vmax)
    ax.axis("off")
    fig.tight_layout(pad=0, h_pad=0, w_pad=0)
    function update(frame)
        im1.set_data(imgs[frame%2+1])
        return [im1]
    end 
    ani = animation.FuncAnimation(fig, update, frames=nFrames, interval=interval, blit=true)
    return fig, ani
end