using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)
BHO = BlochHighOrder("000")

hoseq = demo_hoseq();# plot_hoseqd(hoseq);
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral/results_Stitched_B0"; if ispath(dir) == false mkdir(dir) end
TE = 0.0149415; # s
Nx = Ny = 150;
############################################################################################## 
# Simu
############################################################################################## 
# 1. hoseq

# 2. phantom
obj = brain_phantom2D(brain2D(); ss=simtype.ss, location=0.8, B0_type=:real, B0_file=:B0_medianfiltered_r4); info(obj)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

p_Δw = plot_phantom_map(obj, :Δw; darkmode=true)
# savefig(p_Δw, dir*"/quadraticB0map_$(maxOffresonance)_objΔw.svg", width=500,height=500,format="svg")
ref = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                B0_type=:real, B0_file=:B0_medianfiltered_r4, key=:Δw); 
p_Δw_ref = plot_image(ref; title="B0map, [$(minimum(ref)),$(maximum(ref))] Hz", darkmode=true, zmin=minimum(ref))
savefig(p_Δw_ref, dir*"/realB0map_B0map.svg", width=550,height=500,format="svg")


# 3. scanner & sim_params
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
img = reconstruct_2d_image(raw);
p_image = plot_image(img; darkmode=true, title="Sim: $(BHO.name), Δw: [$(minimum(ref)),$(maximum(ref))] Hz")
savefig(p_image, dir*"/realB0map_reconNUFFT$(BHO.name).svg", width=550,height=500,format="svg")


############################################################################################## 
# Recon
############################################################################################## 

Nx, Ny = raw.params["reconSize"][1:2];
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

B0map = brain_phantom2D_reference(brain2D(); ss=simtype.ss, location=0.8,target_fov=(150, 150), target_resolution=(1,1),
                                    B0map=:file,key=:Δw);
p_ref_B0map = plot_image(B0map; title="B0map, [$(minimum(ref)),$(maximum(ref))] Hz", zmin=minimum(B0map))

#######################################################################################
# iterative SignalOp 
#######################################################################################
["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
reg = "L2"
@info "reg: $reg"
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = reg
recParams[:λ] = 1e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"

Op = HighOrderOp(shape, tr_nominal, tr_skope, BHO; Nblocks=9, fieldmap=Matrix(B0map))
# Op = SignalOp(shape, tr_nominal, 1e-3, 1e-3; Nblocks=9, fieldmap=B0map)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp $(BHO.name) with B0map [$(minimum(ref)),$(maximum(ref))] Hz", width=650, height=600)
savefig(p_iter_SignalOp, dir*"/realB0map_reconHighOrderOp$(BHO.name).svg", width=550,height=500,format="svg")


