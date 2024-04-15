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
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral"

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

mat_111 = MAT.matread("$dir/HighOrderOp_Simu_111.mat")  # keys: ["imgs", "imgs_error", "BHO"]
mat_000 = MAT.matread("$dir/HighOrderOp_Simu_000.mat")



imgs_111    = mat_111["imgs"];
imgs_error_111 = mat_111["imgs_error"];
BHO_recos = mat_111["BHO"];
subplot_titles = ["Reco: $t" for t in BHO_recos];
title="HighOrderOp, Simu: 111";

# p_111       = plot_imgs(imgs_111, subplot_titles; title=title, width=1300, height=200)
# savefig(p_000,       dir*"/HighOrderOp_Simu_000.svg", width=1300, height=200,format="svg")

imgs_error_111 = Array{Float32,3}(undef, size(imgs_111));
imgs_normalized_111 = Array{Float32,3}(undef, size(imgs_111));
for idx in eachindex(BHO_recos)
    imgs_error_111[:,:, idx] = ρ2 - normalization(imgs_111[:,:, idx]);
    imgs_normalized_111[:,:, idx] = normalization(imgs_111[:,:, idx]);
end

p_111_normalized = plot_imgs(imgs_normalized_111, subplot_titles; title=title*", normalized to [0,1]", width=1300, height=200)
p_111_error = plot_imgs(imgs_error_111, subplot_titles; title=title*", error map", width=1300, height=200)
savefig(p_111_error,  dir*"/HighOrderOp_Simu_111_ErrorWithPhantom.svg", width=1300, height=200,format="svg")
savefig(p_111_normalized,  dir*"/HighOrderOp_Simu_111_Normalized.svg", width=1300, height=200,format="svg")

error = standardization(ρ') - standardization(normalize(imgs_111[:,:,8]));
error = ρ - imgs_111[:,:,8];
error = standardization(ρ2) - standardization(normalize(imgs_111[:,:,8]));
plot_image(error; title="$(BHO_recos[8]) - ρ", zmin=minimum(error), zmax=maximum(error))

