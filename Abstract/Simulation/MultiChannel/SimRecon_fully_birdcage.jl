using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances

outpath = "$(@__DIR__)/Abstract/Simulation/MultiChannel/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory
############################################################################################## 
# Setup
############################################################################################## 
T = Float64;
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
location = 0.8;                                              
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
Nx = Ny = 150;
# setting the coil sensitivity used in the simulation
csm_type  = :birdcage;      # a simulated birdcage coil-sensitivity
csm_nCoil = 8;              # 9-channel
csm_nRow  = 2;
csm_nCol  = 4;

db0_type  = :quadratic;     
db0_max   = :100.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map


#############################################################################
# 1. Simulation
#############################################################################
# 1. sequence
hoseq_stitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
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
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft, dims=[1,2]); dim=3, nRow=csm_nRow, nCol=csm_nCol)


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
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft, dims=[1,2]); dim=3, nRow=csm_nRow, nCol=csm_nCol)

#############################################################################
# 3. Preparing for reconstruction
#############################################################################
# Coil-Sensitivity Map
coil = csm_Birdcage(217, 181, csm_nCoil, verbose=true);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, csm_nCoil);
for c = 1:csm_nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
csm = mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2])
fig_csm = plt_images(csm; dim=3, nRow=csm_nRow, nCol=csm_nCol)

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (1., 1.); location=location, ss=simtype.ss, db0_type=db0_type, db0_max=db0_max);
plt_B0map(rotl90(B0map))
x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=simtype.ss);
plt_image(rotl90(x_ref))

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
solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;
Δx = Δy = 1e-3;
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
# recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true);
recParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));

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
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (1., 1.); location=location, ss=simtype.ss, db0_type=db0_type, db0_max=db0_max);
B0map = rotl90(B0map);

x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=simtype.ss);
x_ref = rotl90(x_ref);

headmask = brain_phantom2D_reference(phantom, :headmask, (150., 150.), (1., 1.); location=location, ss=simtype.ss);
headmask = rotl90(headmask);

MAT.matwrite("$(outpath)/fully_snr10_$(solver)_$(iter)_$(regularization)_$(λ).mat", Dict("imgs"=>imgs, "labels"=>labels, "csm"=>csm, "signal"=>data, "B0map"=>B0map, "x_ref"=>x_ref, "headmask"=>headmask))
