using MAT, PlotlyJS, Statistics
using ImageTransformations, ImageQualityIndexes, ImageDistances


folder = "woT2B0"   #  "woT2B0", "woB0_wT2"   
dir = "$(@__DIR__)/src/demo/demo_recon/SignalOp_spiral/results/$folder"


mask = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:binary, target_fov=(150, 150), target_resolution=(1,1));
p_ref_mask = plot_image(mask; title="PhantomReference[ $(size(mask)) | 1mm | binarymask ]")
savefig(p_ref_mask,  dir*"/PhantomReference_ss3_location0.8_binary.svg", width=500, height=450,format="svg")

ρ = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
if folder == "woB0_wT2"
    T2map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1)); # plot_image(T2map; title="T2map", width=650, height=600)
    T2map = (T2map.<=46*1e-3) .* Inf .+ T2map; R2map = 1 ./ T2map;
    ρ = normalization(ρ .* exp.(-0.0149415 .* R2map))
    p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ] withT2")
    savefig(p_ref,  dir*"/PhantomReference_ss3_location0.8_rho_withT2.svg", width=500, height=450,format="svg")
elseif folder == "woT2B0"
    p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ]")
    savefig(p_ref,  dir*"/PhantomReference_ss3_location0.8_rho.svg", width=500, height=450,format="svg")
end

mat = MAT.matread("$dir/SignalOp_Simu_000.mat");  # keys: ["img_direct_NUFFTOp","img_iter_NUFFTOp","img_iter_SignalOp_normalized","img_iter_SignalOp","img_iter_HighOrderOp"]
img_direct_NUFFTOp = mat["img_direct_NUFFTOp"];
img_iter_NUFFTOp = mat["img_iter_NUFFTOp"];
img_iter_SignalOp_normalized = mat["img_iter_SignalOp_normalized"];
img_iter_SignalOp = mat["img_iter_SignalOp"];
img_iter_HighOrderOp = mat["img_iter_HighOrderOp"];

imgs = Array{Float32,3}(undef, size(img_direct_NUFFTOp)..., 5);
imgs[:,:, 1] = img_direct_NUFFTOp; imgs[:,:, 2] = img_iter_NUFFTOp; imgs[:,:, 3] = img_iter_SignalOp_normalized; imgs[:,:, 4] = img_iter_SignalOp; imgs[:,:, 5] = img_iter_HighOrderOp;
p_imgs = plot_imgs(imgs, ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]; title="", width=1100, height=250)

# error map & imgs normalized to [0,1]
imgs_error = Array{Float32,3}(undef, size(imgs));
imgs_normalized = Array{Float32,3}(undef, size(imgs));
imgs_error_mask = Array{Float32,3}(undef, size(imgs));
imgs_normalized_mask = Array{Float32,3}(undef, size(imgs));
for idx = 1:5
    imgs_error[:,:, idx] = ρ - normalization(imgs[:,:, idx]);
    imgs_normalized[:,:, idx] = normalization(imgs[:,:, idx]);
    imgs_error_mask[:,:, idx]  = imgs_error[:,:, idx] .* mask;
    imgs_normalized_mask[:,:, idx] = imgs_normalized[:,:, idx] .* mask;
end


title="Simu: 000, $folder";
subplot_titles = ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]
width=1000; height=210;

# plot_imgs: error map
mse_values = Vector{Float32}(undef, size(imgs)[end])
ssim_values = Vector{Float32}(undef, size(imgs)[end])
for idx = 1:5
    mse_values[idx] = mse(normalization(imgs[:,:, idx]), ρ)
    ssim_values[idx] = assess_ssim(normalization(imgs[:,:, idx]), ρ)
end
annotations = []
for idx in eachindex(subplot_titles)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_error = plot_imgs(imgs_error, subplot_titles; title=title*" | error map", 
                    width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_error,  dir*"/SignalOp_Simu_000_ErrorWithPhantom.svg",format="svg", width=width+100, height=height+40+40)
p_error_abs = plot_imgs(abs.(imgs_error), subplot_titles; title=title*" | error map", 
                    width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_error_abs,  dir*"/SignalOp_Simu_000_ErrorWithPhantom_abs.svg",format="svg", width=width+100, height=height+40+40)


# plot_imgs: error map masked
mse_values = Vector{Float32}(undef, size(imgs)[end])
ssim_values = Vector{Float32}(undef, size(imgs)[end])
for idx = 1:5
    mse_values[idx] = mse(normalization(imgs[:,:, idx])[mask.==1], ρ[mask.==1])
    ssim_values[idx] = assess_ssim(normalization(imgs[:,:, idx]).*mask, ρ.*mask)
end
annotations = []
for idx in eachindex(subplot_titles)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_error_mask = plot_imgs(imgs_error_mask, subplot_titles; title=title*" | error map masked", 
                        width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_error_mask,  dir*"/SignalOp_Simu_000_ErrorWithPhantom_masked.svg",format="svg", width=width+100, height=height+40+40)
p_error_mask_abs = plot_imgs(abs.(imgs_error_mask), subplot_titles; title=title*" | error map masked", 
                        width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_error_mask_abs,  dir*"/SignalOp_Simu_000_ErrorWithPhantom_masked_abs.svg",format="svg", width=width+100, height=height+40+40)



# plot_imgs: imgs normalized to [0,1]
p_normalized = plot_imgs(imgs_normalized, subplot_titles; title=title*" | normalized to [0,1]", width=width, height=height)
p_normalized_mask = plot_imgs(imgs_normalized_mask, subplot_titles; title=title*" | normalized to [0,1], masked", width=width, height=height)
savefig(p_normalized,  dir*"/SignalOp_Simu_000_Normalized.svg",format="svg", width=width+100, height=height+40)
savefig(p_normalized_mask,  dir*"/SignalOp_Simu_000_Normalized_masked.svg",format="svg", width=width+100, height=height+40)


# difference map: SignalOp - SignalOp_NormKspace
diff = img_iter_SignalOp - img_iter_SignalOp_normalized
p_diff = plot_image(diff; title="SignalOp - SignalOp NormKspace", width=660, height=600, zmin=minimum(diff), zmax=maximum(diff))
savefig(p_diff,  dir*"/SignalOp-SignalOp_NormKspace.svg",format="svg", width=660, height=600)