using KomaHighOrder, MRIReco, MAT
import RegularizedLeastSquares: SolverInfo

#########################################################################################
# 1. Load sequence and dynamic field data
#     a. Load sequence file (*.seq format from pulseq)
#     b. Load dynamic field data (*.mat format)
#     c. Create a HO_Sequence object as defined in KomaHighOrder.jl by combining the 
#        sequence and dynamic field data.
#########################################################################################

seq_file = "$(@__DIR__)/demo/MultiChannel/xw_sp2d_7T-1mm-200-r3-noSync-fa90.seq"   # *.seq file is the pulseq's sequence file
# under-sampling factor: R = 4
# inplane resolution: 1 mm x 1 mm
# FOV: 200 mm x 200 mm
# readout duration: 29 ms

dfc_file = "$(@__DIR__)/demo/MultiChannel/xw_sp2d_7T-1mm-200-r3.mat"   # *.mat file contains the dynamic field data from both stitching method and the standard method.
# The dynamic field data is stored in the *.mat file with the following keys:
#= 
"dt"                     [s], time interval between two time points.
"delayStandard"          [s], inital delay time in the field monitoring.
"delayStitched"          [s], inital delay time in the field monitoring.
"nSampleAllSegStandard"  number of time points of the standard method.
"nSampleAllSegStitched"  number of time points of the stitched method. A vector, one element for each segment.
"bfieldStandard"         [mT], fields of dynamic measurment (up to 2nd order) of the standard method.
"bfieldStitched"         [mT], fields of dynamic measurment (up to 2nd order) of the stitched method.
"ksphaStandard"          [rad], coefficients of dynamic measurment (up to 2nd order) of the standard method.
"ksphaStitched"          [rad], coefficients of dynamic measurment (up to 2nd order) of the stitched method.
=#

print("$(@__DIR__)")
# 1. read the *.seq file and reverse the sign of the gradient (axis x)
seq = read_seq(seq_file)[3:end-3]; # read_seq function from KomaMRI.jl, return a struct of Sequence
seq.GR[1,:] = -seq.GR[1,:];        # reverse the sign of the gradient (axis x)
plot_seq(seq)                      # plot_seq function from KomaMRI.jl, plot the Sequence


grad = MAT.matread(dfc_file);
Δt = grad["dt"];
nGradSample, nTerm = size(grad["bfieldStitched"])              
dfcStitched = grad["bfieldStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
dfcStandard = grad["bfieldStandard"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
ntStitched  = grad["nSampleAllSegStitched"]
ntStandard  = grad["nSampleAllSegStandard"]
t = Δt * (nGradSample-1);
GR_dfcStitched = reshape([KomaMRIBase.Grad(dfcStitched[idx,:], t, Δt/2, Δt/2, 0) for idx=1:9], :, 1);
GR_dfcStandard = reshape([KomaMRIBase.Grad(dfcStandard[idx,:], t, Δt/2, Δt/2, 0) for idx=1:9], :, 1);


hoseqStitched = HO_Sequence(seq);                       # hoseq, defined in KomaHighOrder.jl
hoseqStandard = HO_Sequence(seq);                       # hoseq, defined in KomaHighOrder.jl
hoseqStitched.GR_dfc[2:4, :] = hoseqStitched.SEQ.GR;    # copy the 1st-order nominal gradient data from the seq object to the hoseq object
hoseqStandard.GR_dfc[2:4, :] = hoseqStandard.SEQ.GR;    # copy the 1st-order nominal gradient data from the seq object to the hoseq object
hoseqStitched.GR_dfc[:,5] = GR_dfcStitched;             # "5" is the index of the readout block in the spiral sequence
hoseqStandard.GR_dfc[:,5] = GR_dfcStandard;             # "5" is the index of the readout block in the spiral sequence
plot_seq(hoseqStitched)
plot_seq(hoseqStandard)


#########################################################################################
# 2. Simulation with dynamic field changes and ΔB0 map (quadratic-shaped)
#     a. setting the coil-sensitivity used in simulation
#     b. we use dynamic fields from field monitoring (stitching)
#     c. recon the signal with NuFFT immeadiately after simulation
#########################################################################################
# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 1        # set phantom down-sample factor to 5
BHO = BlochHighOrder("111", true, true)                          # ["111"] turn on all order terms of dynamic field change. turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D", x=0.2, y=0.2, z=0.2) # setting for Phantom: decide which phantom file to use, loading phantom from src/phantom/mat folder
maxOffresonance = 100.                                            # set the maximum off-resonance frequency in Hz for quadratic B0 map

# setting the coil sensitivity used in the simulation
csmtype= :birdcage;          # a simulated birdcage coil-sensitivity
nCoil=9;                     # 8-channel

# 1. sequence
plot_seq(hoseqStitched)
Nx=Ny=200;     # matrix size for recon

# 2. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=0.65, B0type=:quadratic, maxOffresonance=maxOffresonance, csmtype=csmtype, nCoil=nCoil)
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # set Δw to 0 if B0=false
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation
# plot_phantom_map(obj, :ρ)
# plot_phantom_map_csm(obj, :mag; coil_idx=1)  # [:pha, :real, :imag] plot the magnitude of coil-sensitivity in channel 1.

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params();
sim_params["sim_method"]  = BHO;      # using "BlochHighOrder" for simulation with high-order terms
sim_params["return_type"] = "mat";    # setting with "mat", return the signal data for all channel
sim_params["precision"]   = "f64";
sim_params["Nblocks"] = 1000;

# 4. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);          
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=4, height=4)
fig_cha = plt_images(permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]),width=4, height=4)

# 5. Adding noise to signal data
snr = 15;
data = signal[:,:,1];
nSample, nCha = size(data);

signalAmpl = sum(abs.(data), dims=1)/ nSample;
data = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));

# show the effect of noise
raw = signal_to_raw_data(reshape(data, (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=4, height=4)
fig_cha = plt_images(permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]),width=4, height=4)

#########################################################################################
# 3. Reconstruction with HighOrderOp
#     a. obtain the ΔB₀ map
#     b. construct the encoding operator (HighOrderOp), defined in our KomaHighOrder.jl
#     c. reconstruct the signal with the encoding operator
#########################################################################################
acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

# Coil-Sensitivity Map
coil = csm_Birdcage(217, 181, nCoil, verbose=true);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
csm = permutedims(mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]), [3, 1, 2]);
fig_csm = plt_images(csm,width=4, height=4)

# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
B0map = brain_phantom2D_reference(phantom; ss=ss, location=0.65, target_fov=(200, 200), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
fig_b0map = plt_image(rotl90(B0map))

# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom; ss=ss, location=0.65, key=:ρ, target_fov=(200, 200), target_resolution=(1,1));
fig_ref = plt_image(rotl90(x_ref))

# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
_, K_nominal_adc, _, K_dfc_adc_stitched = get_kspace(hoseqStitched; Δt=1);
_,             _, _, K_dfc_adc_standard = get_kspace(hoseqStandard; Δt=1);

times = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);

# create the trajectories for the nominal and the stitching measurement
tr_nominal      = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_ksphaStitched = Trajectory(K_dfc_adc_stitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_ksphaStandard = Trajectory(K_dfc_adc_standard'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);


########################################################################
# SENSE
########################################################################
solver = "admm"; reg = "TV"; iter = 50; λ = 1e-4;
solver = "cgnr"; reg = "L2"; iter = 50; λ = 1e-3;

T = Float64;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
@time rec = reconstruction(acqData, recParams).data[:,:];
fig = plt_image(rotl90(abs.(rec)); vminp=0, vmaxp=99.9)  


#############################################################################
# HighOrderOp, the extended signal model for high-order terms
#############################################################################
solver = "admm"; reg = "TV"; iter = 50; λ = 1e-4;
solver = "cgnr"; reg = "L2"; iter = 30; λ = 1e-3;

recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
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

# HOOp = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=5, fieldmap=Matrix(B0map), grid=1, verbose=true);  # 
# Op = DiagOp(HOOp, numChan) ∘ S 
# recParams[:encodingOps] = reshape([Op], 1,1);
# @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
# plt_image(rotl90(rec))




Nblocks=5;
##### with ΔB₀
# 1. nominal trajectory, BlochHighOrder("000")
Op1 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 2. stitched trajectory, BlochHighOrder("110")
Op2 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 3. stitched trajectory, BlochHighOrder("111")
Op3 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 4. standard trajectory, BlochHighOrder("111")
Op4 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStandard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);


##### w/o ΔB₀
# 5. nominal trajectory, BlochHighOrder("000")
Op5 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 6. stitched trajectory, BlochHighOrder("110")
Op6 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 7. stitched trajectory, BlochHighOrder("111")
Op7 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 8. standard trajectory, BlochHighOrder("111")
Op8 = HighOrderOp((Nx, Ny), tr_nominal, tr_ksphaStandard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);

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
