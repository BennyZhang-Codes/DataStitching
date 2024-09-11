using KomaHighOrder

dir = "Figures/Fig5/out"; if ispath(dir) == false mkpath(dir) end

simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use

maxOffresonance = 200.        
Nx = Ny = 150;

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
B0map = rotl90(B0map)
# fig_b0map = plt_image(B0map, title="B0map [-$maxOffresonance, $maxOffresonance] Hz")
x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
x_ref = rotl90(x_ref)

headmask = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:headmask , target_fov=(150, 150), target_resolution=(1,1));
headmask = rotl90(headmask)

x, y = size(headmask);
headcountour = zeros(x, y);
α = zeros(x, y);
for i = 1:x, j = 1:x
    if headmask[i,j] == 1 
        if i == 1 || i == x || j == 1 || j == y
            headcountour[i,j] = 1
            α[i,j] = 1
        elseif (headmask[i-1,j] == 0 || headmask[i+1,j] == 0 || headmask[i,j-1] == 0 || headmask[i,j+1] == 0)
            headcountour[i,j] = 1
            α[i,j] = 1
        end
    end
end

img_headcountour = cat(headcountour, headcountour, headcountour, α; dims=3);  # with alpha channel
img_headref = cat(x_ref, x_ref, x_ref, headmask; dims=3);  # with alpha channel

fig, axs = plt.subplots(1, 4, figsize=(10,5))
for ax in axs
    ax.axis("off")
end
axs[1].imshow(x_ref, cmap="gray")
axs[3].imshow(headmask, cmap="gray")
axs[2].imshow(B0map, cmap="jet")
axs[2].imshow(img_headcountour, cmap="gray")
axs[4].imshow(img_headcountour, cmap="gray")
fig.tight_layout(pad=0, w_pad=0, h_pad=0)
fig.savefig("$(dir)/Fig5_groundtruth.png", dpi=300, bbox_inches="tight")



figure_width       = 5/2.54
figure_height      = 4/2.54


fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(figure_width, figure_height), facecolor=color_facecoler)
for ax in axs
    ax.axis("off")
end

cbs = []
for ax in axs
    ai = ax.imshow(B0map, cmap="jet", interpolation="gaussian") # "bilinear", "spline36", "gaussian"
    divider = mpl_axes_grid1.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.02)
    cb = fig.colorbar(ai, fraction=0.1, aspect=20, shrink=1, cax=cax)
    cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,
                    length=ticklength, width=linewidth, pad=pad_labeltick)
    cb.outline.set_visible(false)
    cb.update_ticks()
    cax.set_title("[Hz]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)
    push!(cbs, cb)
end

axs[1].imshow(img_headref, cmap="gray", alpha=0.1)
axs[2].imshow(img_headcountour, cmap="gray", alpha=1)
fig.tight_layout(pad=0, w_pad=0, h_pad=0)
fig.savefig("$(dir)/Fig5_groundtruth_dB0.png", dpi=300, bbox_inches="tight")





figure_width       = 8/2.54
figure_height      = 4/2.54

ρ_values = sort(unique(x_ref))
cm_gray = plt.cm.get_cmap("gray")
cm_gray_ref = cm_gray.from_list("gray_ref",  cm_gray(ρ_values), length(ρ_values))
norm = mcolors.BoundaryNorm([-1; ρ_values] , cm_gray_ref.N)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), facecolor=color_facecoler)
ax.axis("off")

ai = ax.imshow(x_ref, cmap="gray") # "bilinear", "spline36", "gaussian"
divider = mpl_axes_grid1.make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.02)
cax.set_title("[a.u.]", fontsize=fontsize_ticklabel, color=color_label, pad=pad_label)

cb = fig.colorbar(plt.cm.ScalarMappable(cmap=cm_gray_ref, norm=norm), cax=cax, drawedges=true)
cb.minorticks_off()
cb.ax.tick_params(color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel,
                length=ticklength, width=linewidth, pad=pad_labeltick)
cb.set_ticks(ρ_values - (ρ_values - [-1; ρ_values][1:end-1]) / 2)
cb.set_ticklabels([string(round(v, digits=2)) for v in ρ_values])
cb.ax.tick_params(size=0)
cb.outline.set_edgecolor(color_label)
cb.outline.set_linewidth(linewidth)
cb.dividers.set_color(color_label)
cb.dividers.set_linewidth(linewidth)
# cb.update_ticks()

fig.tight_layout(pad=0, w_pad=0, h_pad=0)
fig.savefig("$(dir)/Fig5_groundtruth_rho.png", dpi=300, bbox_inches="tight")
