# julia -t 4
using CUDA
device!(5) 
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
["admm", "cgnr", "fista", "optista", "pogm", "splitBregman"]
for solver in ["pogm",]
    regularizations = solver == "cgnr" ? ["L2"] : ["L2", "L1", "L21", "TV", "Positive", "Proj"]
    for regularization in regularizations
        i = 10
        mses  = zeros(Float32, i, i)
        ssims = zeros(Float32, i, i)
        λs         = [1/10.0^(idx-1) for idx in range(1, i)]
        iterations = [(idx)*10 for idx in range(1, i)]
        λs_str         = [string(1/10.0^(idx-1)) for idx in range(1, i)]
        iterations_str = [string((idx)*10) for idx in range(1, i)]
        for idx_λ in collect(range(1,i))
            λ = λs[idx_λ]
            for idx_iter in collect(range(1,i))
                iteration = iterations[idx_iter]
                recParams = Dict{Symbol,Any}()
                recParams[:reconSize] = (Nx, Ny)  # 150, 150
                recParams[:densityWeighting] = true
                recParams[:reco] = "standard"
                recParams[:regularization] = regularization
                recParams[:λ] = λ
                recParams[:iterations] = iteration
                recParams[:solver] = solver
                BHO_reco = "000"

                imgs = Array{ComplexF32,2}(undef, Nx, Ny);
                @info "solver: $(solver), regularization: $(regularization), λ: $(λ), iteration: $(iteration)"
                Op = HighOrderOp(shape, tr_nominal, tr_dfc, BlochHighOrder(BHO_reco);  Nblocks=9)
                recParams[:encodingOps] = reshape([Op], 1,1)
                @time rec = reconstruction(acqData, recParams);
                img = abs.(rec.data[:,:])

                mses[ idx_λ, idx_iter] = mse(img, ρ)
                ssims[idx_λ, idx_iter] = assess_ssim(img, ρ)
            end
        end
        mse_min , mse_min_idx  = findmin(mses)
        ssim_max, ssim_max_idx = findmax(ssims)
        title_mse  = "MSE min: $(mse_min), λ: $(1/10.0^(mse_min_idx[1]-1)), iterations: $((mse_min_idx[2])*10)"
        title_ssim = "SSIM max: $(ssim_max), λ: $(1/10.0^(ssim_max_idx[1]-1)), iterations: $((ssim_max_idx[2])*10)"
        p_mse  = plot(heatmap(x=iterations_str, y=λs_str, z=mses), Layout(title=title_mse , height=500, width=500, yaxis=attr(title="λ"), xaxis=attr(title="iterations")))
        p_ssim = plot(heatmap(x=iterations_str, y=λs_str, z=ssims), Layout(title=title_ssim, height=500, width=500, yaxis=attr(title="λ"), xaxis=attr(title="iterations")))
        
        savefig(p_mse,  dir*"/$(solver)_$(regularization)_mse.svg" , width=500, height=500,format="svg");
        savefig(p_ssim, dir*"/$(solver)_$(regularization)_ssim.svg", width=500, height=500,format="svg");
    end
end


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

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-5
recParams[:iterations] = 100
recParams[:solver] = "cgnr"
BHO_reco = "000"

imgs = Array{ComplexF32,2}(undef, Nx, Ny);
Op = HighOrderOp(shape, tr_nominal, tr_dfc, BlochHighOrder(BHO_reco);  Nblocks=9)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
img = abs.(rec.data[:,:])

println("MSE: $(mse(img, ρ))")
println("SSIM: $(assess_ssim(img, ρ))")
plot_image(img)
