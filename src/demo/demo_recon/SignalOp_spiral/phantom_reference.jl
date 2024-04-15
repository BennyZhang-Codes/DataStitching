using MAT, PlotlyJS, Statistics,ImageTransformations
function get_center_range(x::Int64, x_range::Int64)
    center = Int64(floor(x/2))
    return center - Int64(ceil(x_range/2))+1 : center + 1 + Int64(ceil(x_range/2))-1
end

# Normalization
function normalization(img::AbstractArray{T,2}) where T<:Real
    return (img.- minimum(img))./ (maximum(img) .- minimum(img))
end
function standardization(img::AbstractArray{T,2}) where T<:Real
    return (img.- mean(img))./ std(img)
end
dir = "$(@__DIR__)/src/demo/demo_recon/SignalOp_spiral"

ρ1 = demo_obj(; axis="axial", ss=1, location=0.8)';
M, N = size(ρ1)
ρ1 = ρ1[get_center_range(M, 750), get_center_range(N, 750)]; # 1mm 
plot_image(ρ1; title="0.2mm iso, ss=1")
ρ1 = imresize(ρ1, (150, 150));
plot_image(ρ1; title="1mm iso, ss=1")

ρ2 = demo_obj(; axis="axial", ss=5, location=0.8)[get_center_range(217, 150), get_center_range(181, 150)]';
plot_image(ρ2; title="1mm iso, ss=5")

plot_image(ρ1 - ρ2; title="difference: ρ1 - ρ2", zmin=minimum(ρ1 - ρ2), zmax=maximum(ρ1 - ρ2))


savefig(plot_image(ρ'; title="Phantom reference"),  dir*"/HighOrderOp_Simu_PhantomReference.svg", width=400, height=350,format="svg")

mat = MAT.matread("$dir/SignalOp_Simu_000.mat");  # keys: ["img_direct_NUFFTOp","img_iter_NUFFTOp","img_iter_SignalOp_normalized","img_iter_SignalOp","img_iter_HighOrderOp"]
img_direct_NUFFTOp = mat["img_direct_NUFFTOp"];
img_iter_NUFFTOp = mat["img_iter_NUFFTOp"];
img_iter_SignalOp_normalized = mat["img_iter_SignalOp_normalized"];
img_iter_SignalOp = mat["img_iter_SignalOp"];
img_iter_HighOrderOp = mat["img_iter_HighOrderOp"];

imgs = Array{Float32,3}(undef, size(img_direct_NUFFTOp)..., 5);
imgs[:,:, 1] = img_direct_NUFFTOp; imgs[:,:, 2] = img_iter_NUFFTOp; imgs[:,:, 3] = img_iter_SignalOp_normalized; imgs[:,:, 4] = img_iter_SignalOp; imgs[:,:, 5] = img_iter_HighOrderOp
p_imgs = plot_imgs(imgs, ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]; title="", width=1100, height=250)


title="Simu: 000";


imgs_error = Array{Float32,3}(undef, size(imgs));
imgs_normalized = Array{Float32,3}(undef, size(imgs));
for idx = 1:5
    imgs_error[:,:, idx] = ρ2 - normalization(imgs[:,:, idx]);
    imgs_normalized[:,:, idx] = normalization(imgs[:,:, idx]);
end

subplot_titles = ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]
width=1100; height=250;
p_normalized = plot_imgs(imgs_normalized, subplot_titles; title=title*", normalized to [0,1]", width=width, height=height)
p_error = plot_imgs(imgs_error, subplot_titles; title=title*", error map", width=width, height=height)
savefig(p_error,  dir*"/SignalOp_Simu_000_ErrorWithPhantom.svg",format="svg", width=width, height=height)
savefig(p_normalized,  dir*"/SignalOp_Simu_000_Normalized.svg",format="svg", width=width, height=height)


diff = img_iter_SignalOp - img_iter_SignalOp_normalized
p_diff = plot_image(diff; title="SignalOp - SignalOp NormKspace", width=660, height=600, zmin=minimum(diff), zmax=maximum(diff))
savefig(p_diff,  dir*"/SignalOp-SignalOp_NormKspace.svg",format="svg", width=660, height=600)