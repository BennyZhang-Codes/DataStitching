using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter

import ImageTransformations: imresize
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=1)
inplane_voxelsize = 0.05
BHO = BlochHighOrder("111", true, true)
key = :fatsatwo  # :fatsatw, :fatsatwo

dir = "$(@__DIR__)/B0/results/B0_dephassing/$(BHO.name)_$(inplane_voxelsize)_$(String(key))_admm_L1_0.001_50"; if ispath(dir) == false mkpath(dir) end
maxOffresonance = 400.
Nx = Ny = 150;

for maxOffresonance in [100., 200., 300., 400.]
############################################################################################## 
# Simu
############################################################################################## 
# 1. hoseq
hoseq_base = demo_hoseq();# plot_hoseqd(hoseq);
hoseq_dict = Dict(:fatsatw =>hoseq_base, :fatsatwo=>hoseq_base[4:end])
hoseq = hoseq_dict[key];
# 2. phantom
obj = brain_hophantom2D(BrainPhantom(prefix="brain3D724", x=inplane_voxelsize, y=inplane_voxelsize, z=0.2); ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation
# 3. B0map
B01 = quadraticFieldmap(217, 181, maxOffresonance)[:,:,1];
c1 = KomaHighOrder.get_center_range(217, Nx);
c2 = KomaHighOrder.get_center_range(181, Ny);
B0map = B01[c1, c2]';
# B0map = brain_phantom2D_reference(BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2); ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
#                                    B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
# plot_image(B0map, darkmode=true, zmin=-maxOffresonance)
# savefig(plot_image(B0map, darkmode=true, zmin=-maxOffresonance), dir*"/quadraticB0map_217_181_center150_150.svg", width=550,height=500,format="svg")

# 4. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"] = "f64"
sim_params["Nblocks"] = 1000
# 5. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img = recon_2d(raw);
p_image = plot_image(img; darkmode=true, title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")
plot_signal(raw)

############################################################################################## 
# Recon
############################################################################################## 
# Nx, Ny = raw.params["reconSize"][1:2];
Nx = Ny = 150;
acqData = AcquisitionData(raw, BlochHighOrder("111"); sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1);
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
tr_skope = Trajectory(K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

#######################################################################################
# iterative HighOrderOp
#######################################################################################
# for λ in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
λ = 1e-3
solver = "admm"
iter = 50
reg = "L1"
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = reg # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
# recParams = merge(defaultRecoParams(), recParams)

Op = HighOrderOp((Nx, Ny), tr_nominal, tr_skope, BHO; Nblocks=9, fieldmap=Matrix(B0map), grid=1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="HighOrderOp $(BHO.name), $(-maxOffresonance) Hz, $(solver)_$(reg)_iter$(iter)_λ$(λ)", width=650, height=600)
savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconHighOrderOp$(BHO.name)_$(solver)_$(reg)_iter$(iter)_λ$(λ).svg", width=550,height=500,format="svg")
end
#######################################################################################
# standard FieldmapNFFTOp
#######################################################################################
nSample, nCoil = size(acqData.kdata[1])
dt = 1e-6
tAQ = (nSample-1) * dt
acqData.traj[1].AQ=tAQ # important for B0 correction
acqData.traj[1].TE= TE   # 0.0149415
acqData.traj[1].times = TE .+ collect(0:dt:tAQ)

acqData.traj[1].times = times;
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:method] = "nfft"
# recParams = merge(defaultRecoParams(), recParams)

recParams[:regularization] = reg # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver

recParams[:correctionMap] = 1im*B0map*2π;
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="FieldmapNFFTOp, $(-maxOffresonance) Hz, $(solver)_$(reg)_iter$(iter)_λ$(λ)", width=650, height=600)
savefig(p_iter_SignalOp, dir*"/quadraticB0map_$(maxOffresonance)_reconFieldmapNFFTOp$(BHO.name)_$(solver)_$(reg)_iter$(iter)_λ$(λ).svg", width=550,height=500,format="svg")
end

