using MAT, PlotlyJS, Statistics,ImageTransformations


folder = "woT2B0"   #  "woT2B0", "woB0_wT2"   
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

imgs_error_111 = Array{Float32,3}(undef, size(imgs_111));
imgs_normalized_111 = Array{Float32,3}(undef, size(imgs_111));
for idx in eachindex(BHO_recos)
    imgs_error_111[:,:, idx] = ρ - normalization(imgs_111[:,:, idx]);
    imgs_normalized_111[:,:, idx] = normalization(imgs_111[:,:, idx]);
end

p_111_normalized = plot_imgs(imgs_normalized_111, subplot_titles; title=title*" | normalized to [0,1]", width=1300, height=200)
p_111_error = plot_imgs(imgs_error_111, subplot_titles; title=title*" | error map", width=1300, height=200)
savefig(p_111_error,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom.svg", width=1300, height=200,format="svg")
savefig(p_111_normalized,  dir*"/HighOrderOp_Simu_111_Normalized.svg", width=1300, height=200,format="svg")

error = standardization(ρ') - standardization(normalize(imgs_111[:,:,8]));
error = ρ - imgs_111[:,:,8];
error = standardization(ρ2) - standardization(normalize(imgs_111[:,:,8]));
plot_image(error; title="$(BHO_recos[8]) - ρ", zmin=minimum(error), zmax=maximum(error))

