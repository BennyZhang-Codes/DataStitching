using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
csmtype= :birdcage; nCoil=8; nrows=2; ncols=4;
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 100.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;

T = Float64;

dir = "workplace/Abstract/Simulation/MultiChannel/out"; if ispath(dir) == false mkpath(dir) end     # output directory


#############################################################################
# 1. Simulation
#############################################################################
# 1. sequence
hoseq_stitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, csmtype=csmtype, nCoil=nCoil, B0type=:quadratic, maxOffresonance=maxOffresonance); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseq_stitched, sys; sim_params);
data = signal[:,:,1];
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54)
fig_cha = plt_images(permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]),width=8, height=8)


#############################################################################
# 2. Adding noise to signal data
#############################################################################
snr = 10;

data = signal[:,:,1];
nSample, nCha = size(data);
signalAmpl = sum(abs.(data), dims=1)/ nSample;
data = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));
# plt.plot(abs.(data[:, 1]), linewidth=0.5)
raw = signal_to_raw_data(reshape(data, (nSample, nCha, 1)), hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54)

imgs_coil = permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]);
fig_cha = plt_images(imgs_coil,width=8, height=8)

#############################################################################
# 3. Preparing for reconstruction
#############################################################################
# Coil-Sensitivity Map
coil = csm_Birdcage(217, 181, nCoil, verbose=true);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
csm = permutedims(mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]), [3, 1, 2]);
fig_csm = plt_images(csm,width=8, height=8)

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
fig_b0map = plt_image(rotl90(B0map), title="B0map [-$maxOffresonance, $maxOffresonance] Hz")
x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));

acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_dfc_adc_stitched = get_kspace(hoseq_stitched; Δt=1);
_, _, _, K_dfc_adc_standard = get_kspace(hoseq_standard; Δt=1);

times = KomaMRIBase.get_adc_sampling_times(hoseq_stitched.SEQ);

tr_nominal          = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_stitched     = Trajectory(K_dfc_adc_stitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_standard     = Trajectory(K_dfc_adc_standard'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);






#############################################################################
# 4. Running reconstruction
#############################################################################
solver = "admm";
reg = "TV";
λ = 1.e-4;
iter=20;

recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
# recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true);
recParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));

numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, recParams);
senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan);
smaps = senseMaps[:,:,1,:];
S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1)

# Op7 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=10, fieldmap=Matrix(B0map), grid=1, verbose=true);
# Op = DiagOp(Op7, numChan) ∘ S 
# recParams[:encodingOps] = reshape([Op], 1,1);
# @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);

# plt_image(rotl90(rec))


Nblocks=10;
##### w/o ΔB₀
# 1. nominal trajectory, BlochHighOrder("000")
Op1 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 2. stitched trajectory, BlochHighOrder("110")
Op2 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 3. stitched trajectory, BlochHighOrder("111")
Op3 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 4. standard trajectory, BlochHighOrder("111")
Op4 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);

##### with ΔB₀
# 5. nominal trajectory, BlochHighOrder("000")
Op5 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 6. stitched trajectory, BlochHighOrder("110")
Op6 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 7. stitched trajectory, BlochHighOrder("111")
Op7 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 8. standard trajectory, BlochHighOrder("111")
Op8 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
Ops = [Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8];


imgs = Array{T,3}(undef, length(Ops), Nx, Ny);
labels = [ "wB0_nominal",  "wB0_stitched_110",  "wB0_stitched_111",  "wB0_standard_111",
          "woB0_nominal", "woB0_stitched_110", "woB0_stitched_111", "woB0_standard_111",];

for idx in eachindex(Ops)
    Op = DiagOp(Ops[idx], numChan) ∘ S 
    recParams[:encodingOps] = reshape([Op], 1,1);
    @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
    imgs[idx, :, :] = rotl90(rec);
    plt_image(rotl90(rec); title=labels[idx])
end


# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1), B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
B0map = rotl90(B0map);

x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
x_ref = rotl90(x_ref);

headmask = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:headmask , target_fov=(150, 150), target_resolution=(1,1));
headmask = rotl90(headmask);

MAT.matwrite(dir*"/fully_snr10_$(solver)_$(iter)_$(reg)_$(λ).mat", Dict("imgs"=>imgs, "labels"=>labels, "csm"=>csm, "signal"=>data, "B0map"=>B0map, "x_ref"=>x_ref, "headmask"=>headmask))
