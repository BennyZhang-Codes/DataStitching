# julia -t 4
# using CUDA
# device!(1) 

using KomaHighOrder
using MRIReco, MRICoilSensitivities, PlotlyJS, MAT, ImageQualityIndexes, ImageDistances
using ProgressMeter
BHO_simu = "000"
folder = "woB0_wT2"   #  "woT2B0"   

dir = "$(@__DIR__)/src/demo/demo_recon/SignalOp_spiral/results/$folder"; if ispath(dir) == false mkdir(dir) end

raw = demo_raw(BHO_simu; folder=folder)
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

T2map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1));
plot_image(T2map; title="T2map", width=650, height=600)
T2map = (T2map.<=46*1e-3) .* Inf .+ T2map;
R2map = 1 ./ T2map
plot_image(exp.(-0.10 .* R2map))


ρ = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
#######################################################################################
# iterative SignalOp 
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
# recParams[:regularization] = "L1"
# recParams[:λ] = 1e-1
recParams[:iterations] = 1
# recParams[:solver] = "admm"


Op = SignalOp(shape, tr_nominal, 1e-3, 1e-3; Nblocks=9, T2map=T2map.*1)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="iterative SignalOp", width=650, height=600)
img_iter_SignalOp = abs.(rec.data[:,:]);

# wt2 = abs.(rec.data[:,:]);

# error = wt2 - wot2
# plot_image(error; zmin=minimum(error))

mse_value = mse(normalization(img_iter_SignalOp), ρ)
ssim_value = assess_ssim(normalization(img_iter_SignalOp), ρ)
