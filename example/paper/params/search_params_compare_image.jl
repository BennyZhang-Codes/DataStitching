# julia -t 4
using CUDA
device!(0)
CUDA.used_memory()
using KomaHighOrder
using MRIReco, MRICoilSensitivities, PlotlyJS, MAT
using ImageQualityIndexes, ImageDistances
simtype = SimType(B0=false, T2=false, ss=5)
BHO = BlochHighOrder("000")
dfc_method = "Stitched"   # :Stitched or :Standard

sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"] = "f64"

raw = demo_raw(BHO; simtype=simtype)
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData = AcquisitionData(raw, BlochHighOrder("000"); sim_params=sim_params);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq(dfc_method=Symbol(dfc_method))
_, K_nominal_adc, _, K_dfc_adc = get_kspace(hoseq; Δt=1)

tr_dfc = Trajectory(K_dfc_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);
ρ = brain_phantom2D_reference(BrainPhantom(); ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
dir = "$(@__DIR__)/params_search"; if ispath(dir) == false mkdir(dir) end


# ρ = brain_phantom2D_reference(BrainPhantom(); ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
# p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ]")

# println("MSE: ", mse(normalization(img), ρ))
# println("SSIM: ", assess_ssim(normalization(img), ρ))
# admm         ["L2", "L1", "L21", "TV", "Positive", "Proj"]
# cgnr         L2
# fista        ["L2", "L1", "L21", "TV", "Positive", "Proj"]
# optista      ["L2", "L1", "L21", "TV", "Positive", "Proj"] 
# pogm         ["L2", "L1", "L21", "TV", "Positive", "Proj"] 
# splitBregman ["L2", "L1", "L21", "TV", "Positive", "Proj"] 

solvers = Dict(
    1=>Dict("solver"=>"cgnr", "regularization"=>"L2", "λ"=>1e-9, "iterations"=>100),
    2=>Dict("solver"=>"admm", "regularization"=>"TV", "λ"=>1e-5, "iterations"=>20),
    3=>Dict("solver"=>"fista", "regularization"=>"L1", "λ"=>1e-9, "iterations"=>30),
    4=>Dict("solver"=>"fista", "regularization"=>"TV", "λ"=>1e-3, "iterations"=>30),
    5=>Dict("solver"=>"optista", "regularization"=>"TV", "λ"=>1e-3, "iterations"=>20),
    6=>Dict("solver"=>"pogm", "regularization"=>"TV", "λ"=>1e-3, "iterations"=>20),
)

images            = Array{Float64, 3}(undef, Nx, Ny, 6);
images_normalized = Array{Float64, 3}(undef, Nx, Ny, 6);
images_error      = Array{Float64, 3}(undef, Nx, Ny, 6);

for idx = 1:6
    params = solvers[idx]
    solver = params["solver"]
    regularization = params["regularization"]
    λ = params["λ"]
    iterations = params["iterations"]
    @info "Running solver: $(solver) | $(regularization) | $(λ) | $(iterations)"

    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)  # 150, 150
    recParams[:densityWeighting] = true
    recParams[:reco] = "standard"
    recParams[:regularization] = regularization
    recParams[:λ] = λ
    recParams[:iterations] = iterations
    recParams[:solver] = solver
    BHO_reco = "000"
    imgs = Array{ComplexF32,2}(undef, Nx, Ny);
    Op = HighOrderOp(shape, tr_nominal, tr_dfc, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    img = abs.(rec.data[:,:])
    images[:, :, idx] = img

    println("MSE: $(mse(img, ρ))")
    println("SSIM: $(assess_ssim(img, ρ))")
    p_img = plot_image(img; title="$(solver) | $(regularization) | $(iterations) | $(λ)")
    savefig(p_img, dir*"/$(solver)_$(regularization).svg", width=450, height=400,format="svg");
end


p_img = plot_image(ρ; title="Phantom Reference: ρ")
savefig(p_img, dir*"/PhantomReference_ρ.svg", width=450, height=400,format="svg");




for idx = 1:6
    params = solvers[idx]
    solver = params["solver"]
    regularization = params["regularization"]
    λ = params["λ"]
    iterations = params["iterations"]
    @info "Running solver: $(solver) | $(regularization) | $(λ) | $(iterations)"

    img = images[:, :, idx]

    println("MSE: $(mse(normalization(img), ρ))")
    println("SSIM: $(assess_ssim(normalization(img), ρ))")

    images_normalized[:, :, idx] = normalization(img)
    images_error[:, :, idx] = abs.(ρ - normalization(img))
    # p_img            = plot_image(img               ; title="$(solver) | $(regularization) | $(iterations) | $(λ)")
    # p_img_normalized = plot_image(normalization(img); title="$(solver) | $(regularization) | $(iterations) | $(λ), normalized")
    # p_img_error      = plot_image(images_error[:, :, idx]; title="$(solver) | $(regularization) | $(iterations) | $(λ), error")
    # savefig(p_img           , dir*"/$(solver)_$(regularization).svg"           , width=450, height=400,format="svg");
    # savefig(p_img_normalized, dir*"/$(solver)_$(regularization)_normalized.svg", width=450, height=400,format="svg");
    # savefig(p_img_error     , dir*"/$(solver)_$(regularization)_error.svg"     , width=450, height=400,format="svg");
end


mse_values  = Vector{Float32}(undef, 6);
ssim_values = Vector{Float32}(undef, 6);
annotations    = []; 
subplot_titles = [];
for idx = 1:6
    mse_values[idx]  = mse(images_normalized[:,:, idx]        , ρ)
    ssim_values[idx] = assess_ssim(images_normalized[:,:, idx], ρ)
    push!(subplot_titles, "$(solvers[idx]["solver"]) | $(solvers[idx]["regularization"])")
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end

width=1050; height=160;
p_image      = plot_imgs(images           , subplot_titles; title="images"    , width=width, height=height)
p_normalized = plot_imgs(images_normalized, subplot_titles; title="normalized", width=width, height=height)
p_error      = plot_imgs(images_error     , subplot_titles; title="error map" , width=width, height=height, annotations=annotations, margin_bottom=40)

savefig(p_image     , dir*"/images.svg"           , format="svg", width=width+100, height=height+80)
savefig(p_normalized, dir*"/images_normalized.svg", format="svg", width=width+100, height=height+80)
savefig(p_error     , dir*"/images_error.svg"     , format="svg", width=width+100, height=height+80)