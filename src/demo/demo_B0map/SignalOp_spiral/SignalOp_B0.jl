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
dir = "$(@__DIR__)/src/demo/demo_B0map/results"; if ispath(dir) == false mkdir(dir) end

for maxOffresonance in [5., 10., 20., 40., 80., 160., 320., 640.]


############################################################################################## 
# Simu
############################################################################################## 
# 1. hoseq

# 2. phantom
obj = brain_phantom2D(brain2D(); ss=3, location=0.8, B0map=:quadratic, maxOffresonance=maxOffresonance); info(obj)
# obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π  cancel Δw
# obj.T2 .= obj.T2 * Inf;   # cancel T2 relaxiation
p_Δw = plot_phantom_map(obj, :Δw; darkmode=true)
savefig(p_Δw, dir*"/quadraticB0map_$(maxOffresonance)_objΔw.svg", width=500,height=500,format="svg")
ref = brain_phantom2D_reference(brain2D(); ss=3, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                   B0map=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
p_Δw_ref = plot_image(ref; title="quadraticB0map, [-$maxOffresonance,$maxOffresonance] Hz", darkmode=true, zmin=-maxOffresonance)
savefig(p_Δw_ref, dir*"/quadraticB0map_$maxOffresonance.svg", width=550,height=500,format="svg")

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
img = reconstruct_2d_image(raw);
p_image = plot_image(img; darkmode=true, title="Sim: 000, Δw: [-$maxOffresonance,$maxOffresonance] Hz")
savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")


############################################################################################## 
# Recon
############################################################################################## 

Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq()
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1)

t_adc = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ)
# times = t_adc .- minimum(t_adc)
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ)
tr_skope = Trajectory(K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

# T2map = brain_phantom2D_reference(brain2D(); ss=3, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1));
# plot_image(T2map; title="T2map", width=650, height=600)
# T2map = (T2map.<=46*1e-3) .* Inf .+ T2map;
# plot_image(exp.(-0.10 ./ T2map))

B0map = brain_phantom2D_reference(brain2D(); ss=3, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                    B0map=:quadratic,key=:Δw, maxOffresonance=maxOffresonance);
# p_ref_B0map = plot_image(B0map; title="quadraticB0map, [-$maxOffresonance,$maxOffresonance] Hz", zmin=-maxOffresonance)

# ρ = brain_phantom2D_reference(brain2D(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));

# plot_image(ρ)

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
img_iter_SignalOp = abs.(rec.data[:,:]);
savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconSignalOp.svg", width=550,height=500,format="svg")
# wt2 = abs.(rec.data[:,:]);

# error = wt2 - wot2
# plot_image(error; zmin=minimum(error))

# mse_value = mse(normalization(img_iter_SignalOp), ρ)
# ssim_value = assess_ssim(normalization(img_iter_SignalOp), ρ)

end