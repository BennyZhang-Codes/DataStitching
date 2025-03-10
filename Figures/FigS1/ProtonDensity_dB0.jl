using MAT
import Statistics: quantile
using PyPlot, PyCall
matplotlib.rcParams["mathtext.default"]
matplotlib.rc("mathtext", default="regular")
mpl_axes_grid1 = pyimport("mpl_toolkits.axes_grid1")
mcolors = matplotlib.colors

outpath = "$(@__DIR__)/Figures/FigS1/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory



T = Float64;
nX = nY = 500; nZ = 1;  # matrix size for recon
Δx = Δy = 0.3e-3; Δz = 2e-3;
# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 1        # set phantom down-sample factor to 3
res = 0.1
location = 0.8;                                              
phantom = BrainPhantom(prefix="brain3D724", x=0.1, y=0.1, z=0.2) # decide which phantom file to use
db0_type  = :quadratic;     
db0_max   = :150.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map


# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
b0map = brain_phantom2D_reference(phantom, :Δw, (T.(nX*Δx*1e3), T.(nY*Δy*1e3)), (T.(res), T.(res)); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
B0map = rotl90(-b0map);
fig_b0map = plt_B0map(B0map)

# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom, :ρ, (T.(nX*Δx*1e3), T.(nY*Δy*1e3)), (T.(res), T.(res)); location=location, ss=ss);
x_ref = rotl90(x_ref);
fig_ref = plt_image(x_ref)

headmask = brain_phantom2D_reference(phantom, :headmask, (T.(nX*Δx*1e3), T.(nY*Δy*1e3)), (T.(res), T.(res)); location=location, ss=ss);
headmask = rotl90(headmask);
fig_headmask = plt_image(headmask)



x, y = size(headmask);
headcountour = zeros(x, y);
α = zeros(x, y);
for i = 1:x, j = 1:x
    if headmask[i,j] == 1 
        if i == 1 || i == x || j == 1 || j == y
            headcountour[i,j] = 1
            α[i,j] = 1
        elseif (headmask[i-1,j] == 0 || headmask[i-1,j+1] == 0 || headmask[i-1,j-1] == 0 || headmask[i+1,j] == 0 || headmask[i+1,j+1] == 0 || headmask[i+1,j-1] == 0 || headmask[i,j-1] == 0 || headmask[i,j+1] == 0)
            headcountour[i,j] = 1
            α[i,j] = 1
        end
    end
end
for i = 1:3
mask = copy(headcountour);
for i = 2:x-1, j = 2:x-1
    if mask[i,j] == 1 
        # if  i == 1 || i == x || j == 1 || j == y
        # elseif (headmask1[i-1,j] == 0 || headmask1[i+1,j] == 0 || headmask1[i,j-1] == 0 || headmask1[i,j+1] == 0)
        headcountour[i-1,j] = 1; α[i-1,j] = 1
        headcountour[i+1,j] = 1; α[i+1,j] = 1
        headcountour[i-1,j+1] = 1; α[i-1,j+1] = 1
        headcountour[i-1,j-1] = 1; α[i-1,j-1] = 1
        headcountour[i,j+1] = 1; α[i,j+1] = 1
        headcountour[i,j-1] = 1; α[i,j-1] = 1
        headcountour[i+1,j+1] = 1; α[i+1,j+1] = 1
        headcountour[i+1,j-1] = 1; α[i+1,j-1] = 1
    end
end
end
img_headcountour = cat(headcountour, headcountour, headcountour, α; dims=3);  # with alpha channel
img_headref = cat(x_ref, x_ref, x_ref, headmask; dims=3);  # with alpha channel
plt_image(headcountour)


matplotlib.rc("mathtext", default="regular")
# matplotlib.rc("figure", dpi=200)
# matplotlib.rc("font", family="Times New Roman")
matplotlib.rc("font", family="Arial")
matplotlib.rcParams["mathtext.default"]
figure_width       = 5/2.53999863
figure_height      = 5/2.53999863
linewidth          = 0.5
ticklength         = 1.5

fontsize_subfigure = 9
fontsize_label     = 7
fontsize_legend    = 7
fontsize_ticklabel = 6

pad_labeltick      = 2
pad_label          = 2
color_facecolor    = "#ffffff"
color_label        = "#000000"

dpi                = 900;
########################################################################
# Off-resonance map
########################################################################
# fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
# ax = axs[1,1]

fig = plt.figure(figsize=(figure_width, figure_height), facecolor=color_facecolor)
ax    = fig.add_axes([0, 0, 0.7, 0.7])
cax   = fig.add_axes([0.72, 0, 0.05*0.7, 0.7])

ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end

ai = ax.imshow(B0map, cmap="jet")#, interpolation="gaussian") # "bilinear", "spline36", "gaussian"
# ax.imshow(img_headref, cmap="gray", alpha=0.1)
ax.imshow(img_headcountour, cmap="gray", alpha=1)
# ax.set_title("Off-resonance", fontsize=fontsize_label, color=color_label)
# divider = mpl_axes_grid1.make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(ai, cax=cax)
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,length=ticklength, width=linewidth, pad=pad_labeltick)
cb.outline.set_visible(false)
cb.update_ticks()
fig.tight_layout(pad=0)


fig.savefig("$(outpath)/B0map.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.0)
fig.savefig("$(outpath)/B0map.svg", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.0)

########################################################################
# Proton density map
########################################################################
# colorbar setting for the ρ map of the phantom
ρ_values = sort(unique(x_ref))
cm_gray = plt.cm.get_cmap("gray")
cm_gray_ref = cm_gray.from_list("gray_ref",  cm_gray(ρ_values), length(ρ_values))
norm = mcolors.BoundaryNorm([-1; ρ_values] , cm_gray_ref.N)

# fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecolor, squeeze=false)
# ax = axs[1,1]

fig = plt.figure(figsize=(figure_width, figure_height), facecolor=color_facecolor)
ax    = fig.add_axes([0, 0, 0.7, 0.7])
cax   = fig.add_axes([0.72, 0, 0.05*0.7, 0.7])

ax.set_facecolor(color_facecolor)
ax.tick_params(axis="both", bottom=false, top=false, left=false, right=false, labelbottom=false, labeltop=false, labelleft=false, labelright=false)
for spine in ax.spines  # "left", "right", "bottom", "top"
    ax.spines[spine].set_visible(false)
end
ax.imshow(x_ref, cmap="gray") # "bilinear", "spline36", "gaussian"
# divider = mpl_axes_grid1.make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.02)
# ax.set_title("Proton density", fontsize=fontsize_label, color=color_label)
cax.set_title("[a.u.]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm_gray_ref, norm=norm), cax=cax, drawedges=true)
cb.minorticks_off()
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel, length=ticklength, width=linewidth, pad=pad_labeltick)
cb.set_ticks(ρ_values - (ρ_values - [-1; ρ_values][1:end-1]) / 2)
cb.set_ticklabels([string(round(v, digits=2)) for v in ρ_values])
cb.ax.tick_params(size=0)
# cb.outline.set_visible(false)
cb.outline.set_edgecolor(color_label)
cb.outline.set_linewidth(linewidth)
cb.dividers.set_color(color_label)
cb.dividers.set_linewidth(linewidth)
fig.tight_layout(pad=0)

fig.savefig("$(outpath)/ProtonDensity.png", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.0)
fig.savefig("$(outpath)/ProtonDensity.svg", dpi=dpi, transparent=true, bbox_inches="tight", pad_inches=0.0)
