using MAT, PlotlyJS, Statistics
using ImageTransformations, ImageQualityIndexes, ImageDistances, ImageMorphology

simtype = SimType(B0=false, T2=false, ss=5)

skope_method = "Stitched"   # :Stitched or :Standard
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral/results_$skope_method/$(simtype.name)"; if ispath(dir) == false mkdir(dir) end


headmask  = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:headmask , target_fov=(150, 150), target_resolution=(1,1));
brainmask = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:brainmask, target_fov=(150, 150), target_resolution=(1,1));
# brainmask = dilate(brainmask; r=5)
p_ref_headmask  = plot_image(headmask ; title="PhantomReference[ $(size(headmask)) | 1mm | headmask ]");
p_ref_brainmask = plot_image(brainmask; title="PhantomReference[ $(size(brainmask)) | 1mm | brainmask ]");
savefig(p_ref_headmask ,  dir*"/PhantomReference_ss$(simtype.ss)_location0.8_headmask.svg" , width=500, height=450,format="svg");
savefig(p_ref_brainmask,  dir*"/PhantomReference_ss$(simtype.ss)_location0.8_brainmask.svg", width=500, height=450,format="svg");

ρ = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
if simtype.T2 == true
    T2map = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1)); # plot_image(T2map; title="T2map", width=650, height=600)
    T2map = (T2map.<=46*1e-3) .* Inf .+ T2map; R2map = 1 ./ T2map;
    ρ = normalization(ρ .* exp.(-0.0149415 .* R2map))
    p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ] withT2")
elseif simtype.T2 == false
    p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ]")
end
savefig(p_ref,  dir*"/PhantomReference_ss$(simtype.ss)_location0.8_rho.svg", width=500, height=450,format="svg")

mat_111 = MAT.matread("$dir/HighOrderOp_Simu_111.mat");  # keys: ["imgs", "imgs_error", "BHO"]
imgs       = mat_111["imgs"];
imgs_error = mat_111["imgs_error"];
BHO_recos  = mat_111["BHO"];
subplot_titles = ["Reco: $t" for t in BHO_recos];
title = "HighOrderOp, Simu: 111, $(simtype.name)";


# error map & imgs normalized to [0,1]
imgs_headmask            = Array{Float32,3}(undef, size(imgs));
imgs_brainmask           = Array{Float32,3}(undef, size(imgs));
imgs_normalize           = Array{Float32,3}(undef, size(imgs));
imgs_normalize_headmask  = Array{Float32,3}(undef, size(imgs));
imgs_normalize_brainmask = Array{Float32,3}(undef, size(imgs));
imgs_error               = Array{Float32,3}(undef, size(imgs));
imgs_error_headmask      = Array{Float32,3}(undef, size(imgs));
imgs_error_brainmask     = Array{Float32,3}(undef, size(imgs));

for idx in eachindex(BHO_recos)
    imgs_headmask[:,:, idx]            = imgs[:,:, idx] .* headmask;
    imgs_brainmask[:,:, idx]           = imgs[:,:, idx] .* brainmask;
    imgs_normalize[:,:, idx]           = normalization(imgs[:,:, idx]);
    imgs_normalize_headmask[:,:, idx]  = imgs_normalize[:,:, idx] .* headmask;
    imgs_normalize_brainmask[:,:, idx] = imgs_normalize[:,:, idx] .* brainmask;
    imgs_error[:,:, idx]               = ρ - normalization(imgs[:,:, idx]);
    imgs_error_headmask[:,:, idx]      = imgs_error[:,:, idx] .* headmask;
    imgs_error_brainmask[:,:, idx]     = imgs_error[:,:, idx] .* brainmask;
end


mse_values            = Vector{Float32}(undef, length(BHO_recos));
ssim_values           = Vector{Float32}(undef, length(BHO_recos));
mse_values_headmask   = Vector{Float32}(undef, length(BHO_recos));
ssim_values_headmask  = Vector{Float32}(undef, length(BHO_recos));
mse_values_brainmask  = Vector{Float32}(undef, length(BHO_recos));
ssim_values_brainmask = Vector{Float32}(undef, length(BHO_recos));
annotations           = []; 
annotations_headmask  = [];
annotations_brainmask = [];
for idx in eachindex(BHO_recos)
    mse_values[idx]            = mse(imgs_normalize[:,:, idx]               , ρ               )
    mse_values_headmask[idx]   = mse(imgs_normalize[:,:, idx][headmask.==1] , ρ[headmask.==1] )
    mse_values_brainmask[idx]  = mse(imgs_normalize[:,:, idx][brainmask.==1], ρ[brainmask.==1])
    ssim_values[idx]           = assess_ssim(imgs_normalize[:,:, idx]           , ρ             )
    ssim_values_headmask[idx]  = assess_ssim(imgs_normalize[:,:, idx].*headmask , ρ .* headmask )
    ssim_values_brainmask[idx] = assess_ssim(imgs_normalize[:,:, idx].*brainmask, ρ .* brainmask)
    push!(annotations          , attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])"                    ,yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
    push!(annotations_headmask , attr(text="SSIM: $(ssim_values_headmask[idx])<br>MSE: $(mse_values_headmask[idx])"  ,yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
    push!(annotations_brainmask, attr(text="SSIM: $(ssim_values_brainmask[idx])<br>MSE: $(mse_values_brainmask[idx])",yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end

width=1200; height=160;
p_error               = plot_imgs(imgs_error                , subplot_titles; title=title*" | error map"                  , width=width, height=height, annotations=annotations          , margin_bottom=40)
p_error_headmask      = plot_imgs(imgs_error_headmask       , subplot_titles; title=title*" | error map, head mask"       , width=width, height=height, annotations=annotations_headmask , margin_bottom=40)
p_error_brainmask     = plot_imgs(imgs_error_brainmask      , subplot_titles; title=title*" | error map, brain mask"      , width=width, height=height, annotations=annotations_brainmask, margin_bottom=40)
p_error_abs           = plot_imgs(abs.(imgs_error)          , subplot_titles; title=title*" | error map"                  , width=width, height=height, annotations=annotations          , margin_bottom=40)
p_error_headmask_abs  = plot_imgs(abs.(imgs_error_headmask) , subplot_titles; title=title*" | error map, head mask"       , width=width, height=height, annotations=annotations_headmask , margin_bottom=40)
p_error_brainmask_abs = plot_imgs(abs.(imgs_error_brainmask), subplot_titles; title=title*" | error map, brain mask"      , width=width, height=height, annotations=annotations_brainmask, margin_bottom=40)
p_normalize           = plot_imgs(imgs_normalize            , subplot_titles; title=title*" | normalize [0,1]"            , width=width, height=height)
p_normalize_headmask  = plot_imgs(imgs_normalize_headmask   , subplot_titles; title=title*" | normalize [0,1], head mask" , width=width, height=height)
p_normalize_brainmask = plot_imgs(imgs_normalize_brainmask  , subplot_titles; title=title*" | normalize [0,1], brain mask", width=width, height=height)

savefig(p_error              , dir*"/normalize_Simu_111_error.svg"              , format="svg", width=width+100, height=height+80)
savefig(p_error_headmask     , dir*"/normalize_Simu_111_error_headmask.svg"     , format="svg", width=width+100, height=height+80)
savefig(p_error_brainmask    , dir*"/normalize_Simu_111_error_brainmask.svg"    , format="svg", width=width+100, height=height+80)
savefig(p_error_abs          , dir*"/normalize_Simu_111_error_abs.svg"          , format="svg", width=width+100, height=height+80)
savefig(p_error_headmask_abs , dir*"/normalize_Simu_111_error_headmask_abs.svg" , format="svg", width=width+100, height=height+80)
savefig(p_error_brainmask_abs, dir*"/normalize_Simu_111_error_brainmask_abs.svg", format="svg", width=width+100, height=height+80)
savefig(p_normalize          , dir*"/normalize_Simu_111_normalize.svg"          , format="svg", width=width+100, height=height+40)
savefig(p_normalize_headmask , dir*"/normalize_Simu_111_normalize_headmask.svg" , format="svg", width=width+100, height=height+40)
savefig(p_normalize_brainmask, dir*"/normalize_Simu_111_normalize_brainmask.svg", format="svg", width=width+100, height=height+40)

# p_imgs = plot_imgs(imgs, subplot_titles; title=title*" | imgs", width=width, height=height)
# p_imgs_mask = plot_imgs(imgs_mask, subplot_titles; title=title*" | imgs, masked", width=width, height=height)
