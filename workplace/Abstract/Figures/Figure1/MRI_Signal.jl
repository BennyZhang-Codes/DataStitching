using PyPlot
using KomaHighOrder

sh = SphericalHarmonics()

using PyPlot, MAT
using KomaHighOrder
sh = SphericalHarmonics()

path         = "E:/pulseq/20241010_skope_fa90/invivo"
r=4
syn_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin("r$(r)", f)][1]
outpath = "workplace/Abstract/Figures/Figure1/out"; if ispath(outpath) == false mkpath(outpath) end


data       = matread(syn_file)["data"];
nSample, nCha = size(data);


for cha = 1 : nCha
    matplotlib.rc("mathtext", default="regular")
    matplotlib.rc("figure", dpi=200)
    matplotlib.rc("font", family="Times New Roman")
    matplotlib.rcParams["mathtext.default"]
    figure_width       = 3.5/2.54
    figure_height      = 1.8/2.54
    linewidth          = 0.5
    ticklength         = 1.5
    fontsize_legend    = 5
    fontsize_label     = 6
    fontsize_ticklabel = 4
    fontsize_subfigure = 8
    pad_labeltick      = 2
    pad_label          = 2
    color_facecolor    = "#ffffff"
    color_label        = "#000000"
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor)

    ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
    for spine in ax.spines  # "left", "right", "bottom", "top"
        ax.spines[spine].set_color(color_label)
        ax.spines[spine].set_visible(false)
    end
    ax.set_facecolor(color_facecolor)
    ax.plot(abs.(data[:, cha]), linewidth=0.5, color="C$(cha%9)")
    fig.tight_layout(pad=0)

    fig.savefig("$(outpath)/Fig_signal_r$(r)_cha$(cha).png", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
    fig.savefig("$(outpath)/Fig_signal_r$(r)_cha$(cha).svg", dpi=300, bbox_inches="tight", transparent=true, pad_inches=0)
end