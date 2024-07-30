# julia -t 4
# using CUDA
# device!(1) 

using KomaHighOrder
using MRIReco, MRICoilSensitivities, PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter
############################################################################################## 
# Setup
############################################################################################## 
hoseq = demo_hoseq();# plot_hoseqd(hoseq);
maxOffresonance = 5.
dir = "$(@__DIR__)/src/demo/demo_B0map/SignalOp_spiral/results"; if ispath(dir) == false mkdir(dir) end
maxOffresonances = [5., 10., 20., 40., 80., 160., 320., 640.];
TE = 0.0149415; # s
imgs_SignalOp = Array{Float32,3}(undef, Nx, Ny, length(maxOffresonances));
imgs_NUFFT    = Array{Float32,3}(undef, Nx, Ny, length(maxOffresonances));
for maxOffresonance in maxOffresonances


############################################################################################## 
# Simu
############################################################################################## 
# 1. hoseq

# 2. phantom
obj = brain_phantom2D(BrainPhantom(); ss=3, location=0.8, B0map=:quadratic, maxOffresonance=maxOffresonance); info(obj)
# obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π  cancel Δw
# obj.T2 .= obj.T2 * Inf;   # cancel T2 relaxiation
# p_Δw = plot_phantom_map(obj, :Δw; darkmode=true)
# savefig(p_Δw, dir*"/quadraticB0map_$(maxOffresonance)_objΔw.svg", width=500,height=500,format="svg")
ref = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                   B0map=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
# p_Δw_ref = plot_image(ref; title="quadraticB0map, [-$maxOffresonance,$maxOffresonance] Hz", darkmode=true, zmin=-maxOffresonance)
# savefig(p_Δw_ref, dir*"/quadraticB0map_$maxOffresonance.svg", width=550,height=500,format="svg")

# 3. scanner & sim_params
sys = Scanner();
sim_method::BlochHighOrder=BlochHighOrder("000");
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = sim_method;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal);
img = recon_2d(raw);
imgs_NUFFT[:,:,findall(x->x==maxOffresonance, maxOffresonances)] = img
p_image = plot_image(img; darkmode=true, title="Sim: 000, Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")


############################################################################################## 
# Recon
############################################################################################## 

Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq()
_, K_nominal_adc, _, K_dfc_adc = get_kspace(hoseq; Δt=1)

t_adc = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ)
# times = t_adc .- minimum(t_adc)
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ)
tr_dfc = Trajectory(K_dfc_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);


B0map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                    B0map=:quadratic,key=:Δw, maxOffresonance=maxOffresonance);
# p_ref_B0map = plot_image(B0map; title="quadraticB0map, [-$maxOffresonance,$maxOffresonance] Hz", zmin=-maxOffresonance)


# reference with T2 relaxiation
ρ = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
T2map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1));
# plot_image(T2map; title="T2map", width=650, height=600)
T2map = (T2map.<=46*1e-3) .* Inf .+ T2map; R2map = 1 ./ T2map;
ρ = normalization(ρ .* exp.(-0.0149415 .* R2map))
p_ref = plot_image(ρ; title="PhantomReference[ $(size(ρ)) | 1mm | ρ ] withT2")
savefig(p_ref,  dir*"/PhantomReference_ss3_location0.8_rho_withT2.svg", width=500, height=450,format="svg")

#######################################################################################
# iterative SignalOp 
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"

Op = SignalOp(shape, tr_nominal, 1e-3, 1e-3; Nblocks=9, fieldmap=B0map)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="SignalOp with B0map [-$maxOffresonance,$maxOffresonance] Hz, cgnr[λ=1e-2, 20 iterations]", width=650, height=600)
# savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconSignalOp.svg", width=550,height=500,format="svg")

imgs_SignalOp[:,:,findall(x->x==maxOffresonance, maxOffresonances)] = abs.(rec.data[:,:]);
end

# error map & imgs normalized to [0,1]
imgs_SignalOp_error = Array{Float32,3}(undef, size(imgs));
imgs_SignalOp_normalized = Array{Float32,3}(undef, size(imgs));
imgs_NUFFT_normalized = Array{Float32,3}(undef, size(imgs));
for idx in eachindex(BHO_recos)
    imgs_SignalOp_error[:,:, idx] = ρ - normalization(imgs[:,:, idx]);
    imgs_SignalOp_normalized[:,:, idx] = normalization(imgs[:,:, idx]);
    imgs_NUFFT_normalized[:,:, idx] = normalization(imgs_NUFFT[:,:, idx]);
end

#######################################################################################
# save results
#######################################################################################
MAT.matwrite(dir*"/SignalOp_Simu_000_diffB0map.mat", Dict("imgs_NUFFT"=>imgs_NUFFT, "imgs_SignalOp"=>imgs_SignalOp, "maxOffresonances"=>maxOffresonances))

width=1200
height=160
subplot_titles = ["[-$maxOffresonance,$maxOffresonance] Hz" for maxOffresonance in maxOffresonances]
p_imgs_SignalOp  = plot_imgs(abs.(imgs_SignalOp), subplot_titles; title="Simu: 000, with B0, T2 | SignalOp", width=width, height=height)
p_imgs_SignalOp_normalized = plot_imgs(imgs_SignalOp_normalized, subplot_titles; title="Simu: 000, with B0, T2 | SignalOp, normalized to [0,1]", width=width, height=height)
p_imgs_NUFFT = plot_imgs(abs.(imgs_NUFFT), subplot_titles; title="Simu: 000, with B0, T2 | NUFFT", width=width, height=height)
p_imgs_NUFFT_normalized = plot_imgs(abs.(imgs_NUFFT_normalized), subplot_titles; title="Simu: 000, with B0, T2 | NUFFT, normalized to [0,1]", width=width, height=height)

savefig(p_imgs_SignalOp,  dir*"/Simu_000_diffB0map_reconSignalOp.svg",format="svg", width=width+100, height=height+40)
savefig(p_imgs_SignalOp_normalized,  dir*"/Simu_000_diffB0map_reconSignalOp_normalized.svg",format="svg", width=width+100, height=height+40)
savefig(p_imgs_NUFFT,  dir*"/Simu_000_diffB0map_reconNUFFT.svg",format="svg", width=width+100, height=height+40)
savefig(p_imgs_NUFFT_normalized,  dir*"/Simu_000_diffB0map_reconNUFFT_normalized.svg",format="svg", width=width+100, height=height+40)

# plot_imgs: error map
mse_values = Vector{Float32}(undef, length(maxOffresonances))
ssim_values = Vector{Float32}(undef, length(maxOffresonances))
for idx in eachindex(maxOffresonances)
    mse_values[idx] = mse(normalization(imgs[:,:, idx]), ρ)
    ssim_values[idx] = assess_ssim(normalization(imgs[:,:, idx]), ρ)
end
annotations = []
for idx in eachindex(maxOffresonances)
    push!(annotations, attr(text="SSIM: $(ssim_values[idx])<br>MSE: $(mse_values[idx])",
    yanchor="top",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=0,showarrow=false,font=attr(size=14)))
end
p_error = plot_imgs(imgs_SignalOp_error, subplot_titles; title=title*" | error map", 
                        width=width, height=height, annotations=annotations, margin_bottom=40)


