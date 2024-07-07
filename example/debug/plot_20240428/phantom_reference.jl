using KomaHighOrder 
using MAT, PlotlyJS, Statistics
using ImageTransformations, ImageQualityIndexes, ImageDistances

folder = "woT2B0"   #  "woT2B0", "woB0_wT2"   
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral/results/$folder"
out = "$(@__DIR__)/debug_recon/plot_20240428/$folder"; if ispath(out) == false mkdir(out) end

mask = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:binary, target_fov=(150, 150), target_resolution=(1,1));
p_ref_mask = plot_img(mask; title="PhantomReference[ $(size(mask)) | 1mm | binarymask ]", thememode=:light_fontblack, width=500, height=450)
savefig(p_ref_mask,  out*"/PhantomReference_ss3_location0.8_binary.svg", width=500, height=450,format="svg")

ρ = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
if folder == "woB0_wT2"
    T2map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1)); # plot_image(T2map; title="T2map", width=650, height=600)
    T2map = (T2map.<=46*1e-3) .* Inf .+ T2map; R2map = 1 ./ T2map;
    ρ = normalization(ρ .* exp.(-0.0149415 .* R2map))
    p_ref = plot_img(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ] withT2", thememode=:light_fontblack)
    savefig(p_ref,  out*"/PhantomReference_ss3_location0.8_rho_withT2.svg", width=500, height=450,format="svg")
elseif folder == "woT2B0"
    p_ref = plot_img(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ]", thememode=:light_fontblack)
    savefig(p_ref,  out*"/PhantomReference_ss3_location0.8_rho.svg", width=500, height=450,format="svg")
end


mat_111 = MAT.matread("$dir/HighOrderOp_Simu_111.mat")  # keys: ["imgs", "imgs_error", "BHO"]
mat_000 = MAT.matread("$dir/HighOrderOp_Simu_000.mat")


imgs    = mat_111["imgs"];
imgs_error = mat_111["imgs_error"];
BHO_recos = mat_111["BHO"];
subplot_titles = ["Reco: $t" for t in BHO_recos];
title="HighOrderOp, Simu: 111, $folder";


# error map & imgs normalized to [0,1]
imgs_error = Array{Float32,3}(undef, size(imgs));
imgs_normalized = Array{Float32,3}(undef, size(imgs));
imgs_error_mask = Array{Float32,3}(undef, size(imgs));
imgs_normalized_mask = Array{Float32,3}(undef, size(imgs));
for idx in eachindex(BHO_recos)
    imgs_error[:,:, idx] = ρ - normalization(imgs[:,:, idx]);
    imgs_normalized[:,:, idx] = normalization(imgs[:,:, idx]);
    imgs_error_mask[:,:, idx]  = imgs_error[:,:, idx] .* mask;
    imgs_normalized_mask[:,:, idx] = imgs_normalized[:,:, idx] .* mask;
end

width=1200
height=160
# plot_imgs: error map
mse_values = Vector{Float32}(undef, length(BHO_recos))
ssim_values = Vector{Float32}(undef, length(BHO_recos))
for idx in eachindex(BHO_recos)
    mse_values[idx] = mse(normalization(imgs[:,:, idx]), ρ)
    ssim_values[idx] = assess_ssim(normalization(imgs[:,:, idx]), ρ)
end
annotations = []
for idx in eachindex(BHO_recos)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_error = plot_imgs(imgs_error, subplot_titles; title=title*" | error map", 
                        width=width, height=height, annotations=annotations, margin_bottom=40, thememode=:light_fontblack)
p_error_abs = plot_imgs(abs.(imgs_error), subplot_titles; title=title*" | error map", 
                        width=width, height=height, annotations=annotations, margin_bottom=40, thememode=:light_fontblack)
savefig(p_error,  out*"/HighOrderOp_Simu_111_ErrorWithPhantom.svg", width=width+100, height=height+80,format="svg")
savefig(p_error_abs,  out*"/HighOrderOp_Simu_111_ErrorWithPhantom_abs.svg", width=width+100, height=height+80,format="svg")

# plot_imgs: error map masked
mse_values = Vector{Float32}(undef, length(BHO_recos))
ssim_values = Vector{Float32}(undef, length(BHO_recos))
for idx = eachindex(BHO_recos)
    mse_values[idx] = mse(normalization(imgs[:,:, idx])[mask.==1], ρ[mask.==1])
    ssim_values[idx] = assess_ssim(normalization(imgs[:,:, idx]).*mask, ρ.*mask)
end
annotations = []
for idx in eachindex(subplot_titles)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_error_mask = plot_imgs(imgs_error_mask, subplot_titles; title=title*" | error map masked", 
                        width=width, height=height, annotations=annotations, margin_bottom=40, thememode=:light_fontblack)
p_error_mask_abs = plot_imgs(abs.(imgs_error_mask), subplot_titles; title=title*" | error map masked", 
                        width=width, height=height, annotations=annotations, margin_bottom=40, thememode=:light_fontblack)
savefig(p_error_mask,  out*"/HighOrderOp_Simu_111_ErrorWithPhantom_masked.svg",format="svg", width=width+100, height=height+40+40)
savefig(p_error_mask_abs,  out*"/HighOrderOp_Simu_111_ErrorWithPhantom_masked_abs.svg",format="svg", width=width+100, height=height+40+40)


# plot_imgs: imgs normalized to [0,1]
p_normalized = plot_imgs(imgs_normalized, subplot_titles; title=title*" | normalized to [0,1]", width=width, height=height, thememode=:light_fontblack)
p_normalized_mask = plot_imgs(imgs_normalized_mask, subplot_titles; title=title*" | normalized to [0,1], masked", width=width, height=height, thememode=:light_fontblack)
savefig(p_normalized,  out*"/HighOrderOp_Simu_111_Normalized.svg",format="svg", width=width+100, height=height+40)
savefig(p_normalized_mask,  out*"/HighOrderOp_Simu_111_Normalized_masked.svg",format="svg", width=width+100, height=height+40)
