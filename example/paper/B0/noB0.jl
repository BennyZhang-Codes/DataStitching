using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter


import ImageTransformations: imresize
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=false, T2=false, ss=5)
BHO = BlochHighOrder("111")

dir = "$(@__DIR__)/B0"; if ispath(dir) == false mkpath(dir) end
TE = 0.0149415; # s
Nx = Ny = 150;


############################################################################################## 
# Simu
############################################################################################## 
# 1. hoseq
hoseq = demo_hoseq();# plot_hoseqd(hoseq);
# 2. phantom
obj = brain_phantom2D(BrainPhantom(); ss=simtype.ss, location=0.8); info(obj)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation


# 4. scanner & sim_params
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
p_image = plot_image(img; darkmode=true, title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")


############################################################################################## 
# Recon
############################################################################################## 

Nx, Ny = raw.params["reconSize"][1:2];
Nx = Ny = 150
acqData1 = AcquisitionData(raw);
acqData = AcquisitionData(raw, BlochHighOrder("111"); sim_params=sim_params);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq();
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1);

t_adc = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
# times = t_adc .- minimum(t_adc)
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
tr_skope = Trajectory(K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);




#######################################################################################
# iterative SignalOp 
#######################################################################################
["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
reg = "L2"
# for reg in ["L2", "L1", "L21", "TV", "Positive", "Proj"]
@info "reg: $reg"
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = reg
recParams[:λ] = 1e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"

Op = HighOrderOp(shape, tr_nominal, tr_skope, BHO; Nblocks=9, grid=1)
# Op = SignalOp(shape, tr_nominal, 1e-3, 1e-3; Nblocks=9, fieldmap=B0map)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp 111", width=650, height=600)



