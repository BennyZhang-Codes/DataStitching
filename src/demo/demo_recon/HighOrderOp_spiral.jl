# julia -t 4
using CUDA
device!(1) 

using MRIReco, KomaHighOrder, MRICoilSensitivities, ProgressMeter

BHO_simu = "000"
raw = demo_raw(BHO_simu)
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq()
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1)

tr_skope = Trajectory(K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 3
recParams[:solver] = "cgnr"

BHO_reco = "000"
Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=3)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
plot_image(abs.(rec.data[:,:]); title="Simu-$(BHO_simu)_Reco-$(BHO_reco)_CG-HighOrderOp", width=650, height=600)

Op = SignalOp(shape, acqData.traj[1]; Nblocks=3)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
plot_image(abs.(rec.data[:,:]); title="Sim111_CG-SignalOp", width=650, height=600)



#######################################################################################
# direct NUFFT
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
@time rec = reconstruction(acqData, recParams);
plot_image(abs.(rec.data[:,:]); title="direct NUFFT", width=650, height=600)

#######################################################################################
# iterative NUFFTOp
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
@time rec = reconstruction(acqData, recParams);
plot_image(abs.(rec.data[:,:]); title="iterative NUFFTOp", width=650, height=600)
img_iter_NUFFTOp = abs.(rec.data[:,:])

#######################################################################################
# iterative SignalOp
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
Op = SignalOp(shape, acqData.traj[1]; Nblocks=3)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
plot_image(abs.(rec.data[:,:]); title="iterative SignalOp", width=650, height=600)
img_iter_SignalOp = abs.(rec.data[:,:])

plot_image(img_iter_NUFFTOp - img_iter_SignalOp; title="iter_NUFFTOp - iter_SignalOp", width=650, height=600)

#######################################################################################
# iterative HighOrderOp
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
BHO_reco = "000"
Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=3)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
plot_image(abs.(rec.data[:,:]); title="iterative HighOrderOp", width=650, height=600)
img_iter_HighOrderOp = abs.(rec.data[:,:])

plot_image(img_iter_HighOrderOp - img_iter_SignalOp; title="iter_HighOrderOp - iter_SignalOp", width=650, height=600)

#######################################################################################
# HighOrderOp 
# [1] Simu: 111, Reco: 000
# [2] Simu: 111, Reco: 100
# [3] Simu: 111, Reco: 110
# [4] Simu: 111, Reco: 111
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
for BHO_reco in ["000", "100", "110", "111"]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=3)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    plot_image(abs.(rec.data[:,:]); title="Simu-$(BHO_simu)_Reco-$(BHO_reco)_CG-HighOrderOp", width=650, height=600)
end



