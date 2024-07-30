using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using MAT, ImageQualityIndexes, ImageDistances
# using ProgressMeter

import ImageTransformations: imresize
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 200.
Nx = Ny = 150;

dir = "Figures/out"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
hoseq_stitched = demo_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = demo_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img = recon_2d(raw);
p_image = plot_image(img; darkmode=true, title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
plot_image(B0map, darkmode=true, zmin=-maxOffresonance)


acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1);
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
tr_skope   = Trajectory(    K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);



recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = 1e-2
recParams[:iterations] = 30
recParams[:solver] = "cgnr"



# ideally, with BlochHighOrder("000"), we use nominal trajectory and a null B0map
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_skope, BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map).*0, grid=1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp 111 with B0map [-$maxOffresonance,$maxOffresonance] Hz", width=650, height=600)
# savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconHighOrderOp111.svg", width=550,height=500,format="svg")



# include ΔB₀, with BlochHighOrder("000"), we use nominal trajectory
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_skope, BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map).*0, grid=1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp 111 with B0map [-$maxOffresonance,$maxOffresonance] Hz", width=650, height=600)
# savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconHighOrderOp111.svg", width=550,height=500,format="svg")


# include both ΔB₀ and 0th/1st order terms of dynamic field change, with BlochHighOrder("110")
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_skope, BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map).*0, grid=1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp 111 with B0map [-$maxOffresonance,$maxOffresonance] Hz", width=650, height=600)
# savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconHighOrderOp111.svg", width=550,height=500,format="svg")


# include both ΔB₀ and all order terms of dynamic field change, with BlochHighOrder("111")
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_skope, BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map).*0, grid=1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp 111 with B0map [-$maxOffresonance,$maxOffresonance] Hz", width=650, height=600)
# savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconHighOrderOp111.svg", width=550,height=500,format="svg")


