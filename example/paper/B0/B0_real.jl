using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter

import ImageTransformations: imresize
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)
BHO = BlochHighOrder("111", true, true)
key = :fatsatwo  # :fatsatw, :fatsatwo
B0_scale = 100
dir = "$(@__DIR__)/B0/results/B0_real/$(BHO.name)_$(String(key))_admm_L1_0.001_50"; if ispath(dir) == false mkpath(dir) end
Nx = Ny = 150;

for B0_scale in [20, 40, 60, 80, 100]
############################################################################################## 
# Simu
############################################################################################## 
# 1. hoseq
hoseq_base = demo_hoseq();# plot_hoseqd(hoseq);
hoseq_dict = Dict(:fatsatw =>hoseq_base, :fatsatwo=>hoseq_base[4:end])
hoseq = hoseq_dict[key];

# 2. phantom
obj = brain_phantom2D(BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2); ss=simtype.ss, location=0.8, 
                       B0type=:real, B0_file=:B0_medianfiltered_r4); info(obj)
obj.Δw .= simtype.B0 ? obj.Δw*B0_scale/100 : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. B0map
B0map = B0_scale/100*brain_phantom2D_reference(BrainPhantom(); ss=simtype.ss, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                B0type=:real, B0_file=:B0_medianfiltered_r4, key=:Δw); 
p_B0map = plot_image(B0map; title="B0map, scale: $(B0_scale)%, [$(minimum(B0map)),$(maximum(B0map))] Hz", darkmode=true, zmin=minimum(B0map))
savefig(p_B0map, dir*"/realB0map_scale$(B0_scale)percent_B0map.svg", width=550,height=500,format="svg")# 4. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"] = "f64"
# sim_params["Nblocks"] = 10000
# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img = recon_2d(raw);
p_image = plot_image(img; darkmode=true, title="Sim: $(BHO.name), scale: $(B0_scale)%, [$(minimum(B0map)),$(maximum(B0map))] Hz")
savefig(p_image, dir*"/realB0map_scale$(B0_scale)percent_reconNUFFT_$(BHO.name).svg", width=550,height=500,format="svg")

############################################################################################## 
# Recon
############################################################################################## 
Nx = Ny = 150
# Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw, BlochHighOrder("111"); sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1);
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
tr_skope   = Trajectory(    K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

#######################################################################################
# iterative HighOrderOp
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L1"  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = 1e-3
recParams[:iterations] = 50
recParams[:solver] = "admm"

Op = HighOrderOp((Nx, Ny), tr_nominal, tr_skope, BHO; Nblocks=9, fieldmap=Matrix(B0map), grid=1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp $(BHO.name) with B0map, scale: $(B0_scale)%, [$(minimum(B0map)),$(maximum(B0map))] Hz", width=650, height=600)
savefig(p_iter_SignalOp, dir*"/realB0map_scale$(B0_scale)percent_reconHighOrderOp_$(BHO.name).svg", width=550,height=500,format="svg")

end


