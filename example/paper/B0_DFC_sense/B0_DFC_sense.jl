using KomaHighOrder, MRISampling, MRIReco, MRICoilSensitivities, MRISimulation
import KomaHighOrder.MRIBase: rawdata
R = 2
BHO      = BlochHighOrder("000", true, true)
simtype  = SimType(B0=true, T2=false, ss=1)
ph       = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2)
maxdB0   = 400.
csmtype  = :real_32cha
nCoil    = 32; nrows=4; ncols=8;

folder   = "spiral_$(csmtype)_nCoil$(nCoil)"
path     = "$(@__DIR__)/B0_DFC_sense/$folder"; if ispath(path) == false mkpath(path) end

seq = load_seq(seqname="spiral", r=R)[2:end]
hoseq = HO_Sequence(seq)
Nx = Ny = 150
# plot_seq(hoseq)
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

obj = brain_hophantom2D(ph; ss=simtype.ss, location=0.8, csmtype=csmtype, nCoil=nCoil, B0type=:quadratic, maxOffresonance=maxdB0); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

# simulate
signal = simulate(obj, hoseq, sys; sim_params);

raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "spiral_R$(R)_nCoil$(nCoil)"
raw.params["protocolName"] = filename
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)
images = recon_2d(raw, Nx=Nx, Ny=Ny);

# plot_signal(raw)
p_images = plot_imgs_subplots(abs.(images), nrows, ncols; title="$(nCoil) coils: NUFFT recon", height=400, width=800)
p_sos = plot_image(sqrt.(sum(images.^2; dims=3))[:,:,1]; title="$(nCoil) coils: NUFFT recon, SOS", height=400, width=450)
savefig(p_sos   , "$(path)/$(filename)-nufft_SOS.svg", format="svg", height=400, width=450)
savefig(p_images, "$(path)/$(filename)-nufft.svg"    , format="svg", height=400, width=800)



acqData = AcquisitionData(raw); # raw = RawAcquisitionData(mrd);
acqData.traj[1].circular = false;
p_traj = plot_traj2d(acqData.traj[1]; height=400, width=400)
savefig(p_traj, "$(path)/$(filename)-traj.svg", format="svg", height=400, width=400)

shape = (Nx, Ny);
T = Float32;
#############################################################################
# recon with the coil sensitivities as the same used in the simulation
#############################################################################

# 1. CSM
coil = csmtype == :real_32cha ? csm_Real_32cha(217, 181) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
coil = get_center_crop(coil, Nx, Ny);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
p_smap_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(nCoil) coils: Coil Sensitivity (Simulation), SOS")
p_smap_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")
p_smap_pha = plot_imgs_subplots(angle.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")

# 2. B0
B01 = quadraticFieldmap(217, 181, maxdB0)[:,:,1];
c1 = KomaHighOrder.get_center_range(217, Nx);
c2 = KomaHighOrder.get_center_range(181, Ny);
B0map = B01[c1, c2]';

# 3.
acqData = AcquisitionData(raw, BlochHighOrder("111"); sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_dfc_adc = get_kspace(hoseq; Δt=1);
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
tr_dfc   = Trajectory(    K_dfc_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 20
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 1
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));


numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, params);
senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan);
# ft = SignalOp((N, N), tr; Nblocks=3, use_gpu=true)
ft = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc, BHO; Nblocks=9, fieldmap=Matrix(B0map), grid=1);
smaps = senseMaps[:,:,1,:]
S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1)
Op = DiagOp(ft, numChan) ∘ S 
params[:encodingOps] = reshape([Op], 1,1)


img_recon = Array{ComplexF32,2}(undef, Nx, Ny);
@time img_recon = reconstruction(acqData, params).data;
p_img_recon = plot_image(abs.(img_recon[:,:]), title="$(nCoil) coils: Sense recon, R = $(R), sim-sensitivity")

savefig(p_smap_mag , "$(path)/$(filename)_simSmap_mag.svg"  , format="svg", height=400, width=800)
savefig(p_smap_pha , "$(path)/$(filename)_simSmap_pha.svg"  , format="svg", height=400, width=800)
savefig(p_smap_sos , "$(path)/$(filename)_simSmap_sos.svg"  , format="svg", height=400, width=450)
savefig(p_img_recon, "$(path)/$(filename)-simSmap_Sense.svg", format="svg", height=400, width=450)







