using MAT, PlotlyJS, Statistics
using ImageTransformations, ImageQualityIndexes, ImageDistances

folder = "woB0_wT2"   #  "woT2B0", "woB0_wT2"   
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral/results/$folder"

ρ = brain_phantom2D_reference(brain2D(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
p_ref = plot_image(ρ; title="PhantomReference[$(size(ρ)) | 1mm | ρ]")
savefig(p_ref,  dir*"/PhantomReference_ss3_location0.8.svg", width=400, height=350,format="svg")
mat_111 = MAT.matread("$dir/HighOrderOp_Simu_111.mat")  # keys: ["imgs", "imgs_error", "BHO"]
mat_000 = MAT.matread("$dir/HighOrderOp_Simu_000.mat")


imgs_111    = mat_111["imgs"];
imgs_error_111 = mat_111["imgs_error"];
BHO_recos = mat_111["BHO"];
subplot_titles = ["Reco: $t" for t in BHO_recos];
title="HighOrderOp, Simu: 111, $folder";

# metrics 
mse_values = Vector{Float32}(undef, length(BHO_recos))
ssim_values = Vector{Float32}(undef, length(BHO_recos))
for idx in eachindex(BHO_recos)
    mse_values[idx] = mse(normalization(imgs_111[:,:, idx]), ρ)
    ssim_values[idx] = assess_ssim(normalization(imgs_111[:,:, idx]), ρ)
end

# error map & imgs normalized to [0,1]
imgs_error_111 = Array{Float32,3}(undef, size(imgs_111));
imgs_normalized_111 = Array{Float32,3}(undef, size(imgs_111));
for idx in eachindex(BHO_recos)
    imgs_error_111[:,:, idx] = ρ - normalization(imgs_111[:,:, idx]);
    imgs_normalized_111[:,:, idx] = normalization(imgs_111[:,:, idx]);
end

width=1200
height=160

annotations = []
for idx in eachindex(BHO_recos)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_111_error = plot_imgs(imgs_error_111, subplot_titles; title=title*" | error map", 
                        width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_111_error,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom.svg", width=width+100, height=height+80,format="svg")

p_111_normalized = plot_imgs(imgs_normalized_111, subplot_titles; title=title*" | normalized to [0,1]", width=width, height=height)
savefig(p_111_normalized,  dir*"/HighOrderOp_Simu_111_Normalized.svg", width=width+100, height=height+40,format="svg")


