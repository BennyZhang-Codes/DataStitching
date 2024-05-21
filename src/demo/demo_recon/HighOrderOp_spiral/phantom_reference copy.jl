using MAT, PlotlyJS, Statistics
using ImageTransformations, ImageQualityIndexes, ImageDistances, ImageMorphology

simtype = SimType(B0=false, T2=false, ss=5)

skope_method = "Stitched"   # :Stitched or :Standard
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral/results_$skope_method/$(simtype.name)"; if ispath(dir) == false mkdir(dir) end


mask = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:binary, target_fov=(150, 150), target_resolution=(1,1));
# mask = dilate(mask; r=3)
p_ref_mask = plot_image(mask; title="PhantomReference[ $(size(mask)) | 1mm | binarymask ]")
savefig(p_ref_mask,  dir*"/PhantomReference_ss$(simtype.ss)_location0.8_binary.svg", width=500, height=450,format="svg")

ρ = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
if simtype.T2 == true
    T2map = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1)); # plot_image(T2map; title="T2map", width=650, height=600)
    T2map = (T2map.<=46*1e-3) .* Inf .+ T2map; R2map = 1 ./ T2map;
    ρ = ρ .* exp.(-0.0149415 .* R2map)
    p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ] withT2")
elseif simtype.T2 == false
    p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ]")
end
savefig(p_ref,  dir*"/PhantomReference_ss$(simtype.ss)_location0.8_rho.svg", width=500, height=450,format="svg")

mat_111 = MAT.matread("$dir/HighOrderOp_Simu_111.mat")  # keys: ["imgs", "imgs_error", "BHO"]

imgs    = mat_111["imgs"];
imgs_error = mat_111["imgs_error"];
BHO_recos = mat_111["BHO"];
subplot_titles = ["Reco: $t" for t in BHO_recos];
title="HighOrderOp, Simu: 111, $(simtype.name)";

scales = []
mses = []
for s in range(1, 10, 91)
    # print(s, "  ")
    # println(mse(imgs[:,:, 8]*s, ρ))
    append!(scales, s)
    append!(mses, mse(imgs[:,:, 8]*s, ρ))
end
_, min_idx = findmin(mses); scale = scales[min_idx] 
imgs = imgs .* scale

# error map & imgs normalized to [0,1]
imgs_mask       = Array{Float32,3}(undef, size(imgs));
imgs_error      = Array{Float32,3}(undef, size(imgs));
imgs_error_mask = Array{Float32,3}(undef, size(imgs));
for idx in eachindex(BHO_recos)
    imgs_mask[:,:, idx]  = imgs[:,:, idx] .* mask;
    imgs_error[:,:, idx] = ρ - imgs[:,:, idx];
    imgs_error_mask[:,:, idx]  = imgs_error[:,:, idx] .* mask;
end

width=1200; height=160;
# plot_imgs: error map
mse_values       = Vector{Float32}(undef, length(BHO_recos));
ssim_values      = Vector{Float32}(undef, length(BHO_recos));
mse_values_mask  = Vector{Float32}(undef, length(BHO_recos));
ssim_values_mask = Vector{Float32}(undef, length(BHO_recos));
annotations = []; annotations_mask = []
for idx in eachindex(BHO_recos)
    mse_values[idx]       = mse(imgs[:,:, idx], ρ)
    ssim_values[idx]      = assess_ssim(imgs[:,:, idx], ρ)
    mse_values_mask[idx]  = mse(imgs[:,:, idx][mask.==1], ρ[mask.==1])
    ssim_values_mask[idx] = assess_ssim(imgs[:,:, idx].*mask, ρ.*mask)

    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
    push!(annotations_mask, attr(text="SSIM: $(ssim_values_mask[idx])<br>MSE: $(mse_values_mask[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end

p_error = plot_imgs(imgs_error, subplot_titles; title=title*" | error map, scale $scale", 
                        width=width, height=height, annotations=annotations, margin_bottom=40)
p_error_abs = plot_imgs(abs.(imgs_error), subplot_titles; title=title*" | error map, scale $scale", 
                        width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_error,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom.svg", width=width+100, height=height+80,format="svg")
savefig(p_error_abs,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom_abs.svg", width=width+100, height=height+80,format="svg")

# plot_imgs: error map masked
p_error_mask = plot_imgs(imgs_error_mask, subplot_titles; title=title*" | error map masked, scale $scale", 
                        width=width, height=height, annotations=annotations_mask, margin_bottom=40)
p_error_mask_abs = plot_imgs(abs.(imgs_error_mask), subplot_titles; title=title*" | error map masked, scale $scale", 
                        width=width, height=height, annotations=annotations_mask, margin_bottom=40)
savefig(p_error_mask,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom_masked.svg",format="svg", width=width+100, height=height+40+40)
savefig(p_error_mask_abs,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom_masked_abs.svg",format="svg", width=width+100, height=height+40+40)


p_imgs = plot_imgs(imgs, subplot_titles; title=title*" | images", width=width, height=height)
p_imgs_mask = plot_imgs(imgs_mask, subplot_titles; title=title*" | images, masked", width=width, height=height)
savefig(p_imgs,  dir*"/HighOrderOp_Simu_111.svg",format="svg", width=width+100, height=height+40)
savefig(p_imgs_mask,  dir*"/HighOrderOp_Simu_111_masked.svg",format="svg", width=width+100, height=height+40)
