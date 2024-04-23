using MAT, PlotlyJS, Statistics
using ImageTransformations, ImageQualityIndexes, ImageDistances


folder = "woT2B0"   #  "woT2B0", "woB0_wT2"   
dir = "$(@__DIR__)/src/demo/demo_recon/SignalOp_spiral/results/$folder"

ρ = brain_phantom2D_reference(brain2D(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
ρ = Matrix{Float32}(ρ);
p_ref = plot_image(ρ; title="PhantomReference[$(size(ρ)) | 1mm | ρ]")
savefig(p_ref,  dir*"/PhantomReference_ss3_location0.8.svg", width=400, height=350,format="svg")

mat = MAT.matread("$dir/SignalOp_Simu_000.mat");  # keys: ["img_direct_NUFFTOp","img_iter_NUFFTOp","img_iter_SignalOp_normalized","img_iter_SignalOp","img_iter_HighOrderOp"]
img_direct_NUFFTOp = mat["img_direct_NUFFTOp"];
img_iter_NUFFTOp = mat["img_iter_NUFFTOp"];
img_iter_SignalOp_normalized = mat["img_iter_SignalOp_normalized"];
img_iter_SignalOp = mat["img_iter_SignalOp"];
img_iter_HighOrderOp = mat["img_iter_HighOrderOp"];

imgs = Array{Float32,3}(undef, size(img_direct_NUFFTOp)..., 5);
imgs[:,:, 1] = img_direct_NUFFTOp; imgs[:,:, 2] = img_iter_NUFFTOp; imgs[:,:, 3] = img_iter_SignalOp_normalized; imgs[:,:, 4] = img_iter_SignalOp; imgs[:,:, 5] = img_iter_HighOrderOp;
p_imgs = plot_imgs(imgs, ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]; title="", width=1100, height=250)


# metrics 
mse_values = Vector{Float32}(undef, size(imgs)[end])
ssim_values = Vector{Float32}(undef, size(imgs)[end])
for idx = 1:5
    mse_values[idx] = mse(normalization(imgs[:,:, idx]), ρ)
    ssim_values[idx] = assess_ssim(normalization(imgs[:,:, idx]), ρ)
end

# error map & imgs normalized to [0,1]
imgs_error = Array{Float32,3}(undef, size(imgs));
imgs_normalized = Array{Float32,3}(undef, size(imgs));
for idx = 1:5
    imgs_error[:,:, idx] = ρ - normalization(imgs[:,:, idx]);
    imgs_normalized[:,:, idx] = normalization(imgs[:,:, idx]);
end


title="Simu: 000, $folder";
subplot_titles = ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]
width=1000; height=210;
p_normalized = plot_imgs(imgs_normalized, subplot_titles; title=title*" | normalized to [0,1]", width=width, height=height)
savefig(p_normalized,  dir*"/SignalOp_Simu_000_Normalized.svg",format="svg", width=width+100, height=height+40)

annotations = []
for idx in eachindex(subplot_titles)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_error = plot_imgs(imgs_error, subplot_titles; title=title*" | error map", width=width, height=height, annotations=annotations, margin_bottom=40)
savefig(p_error,  dir*"/SignalOp_Simu_000_ErrorWithPhantom.svg",format="svg", width=width+100, height=height+40+40)


# difference map: SignalOp - SignalOp_NormKspace
diff = img_iter_SignalOp - img_iter_SignalOp_normalized
p_diff = plot_image(diff; title="SignalOp - SignalOp NormKspace", width=660, height=600, zmin=minimum(diff), zmax=maximum(diff))
savefig(p_diff,  dir*"/SignalOp-SignalOp_NormKspace.svg",format="svg", width=660, height=600)