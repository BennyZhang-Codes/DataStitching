using HighOrderMRI, MRIReco, MAT
import RegularizedLeastSquares: SolverInfo

#########################################################################################
# 1. Load sequence and dynamic field data
#     a. Load sequence file (*.seq format from pulseq)
#     b. Load dynamic field data (*.mat format)
#     c. Create a HO_Sequence object as defined in HighOrderMRI.jl by combining the 
#        sequence and dynamic field data.
#########################################################################################
path     = joinpath(@__DIR__, "demo/Sim_MultiChannel")
seq_file = "$(path)/7T_1p0_200_r4.seq"   # *.seq file is the pulseq's sequence file
# under-sampling factor: R = 4
# inplane resolution: 1 mm x 1 mm
# FOV: 200 mm x 200 mm
# readout duration: ~29 ms

dfc_file = "$(path)/7T_1p0_200_r4.mat"   # *.mat file contains the nominal and measured trajectories.

dfc_data = MAT.matread(dfc_file);
dt_Stitched   = dfc_data["dt_Stitched"];
dt_Standard   = dfc_data["dt_Standard"];
dt_Nominal    = dfc_data["dt_Nominal"];

ksphaStitched = -dfc_data["ksphaStitched"];   # coefficients of dynamic measurment (up to 3rd order) of the data stitching method
ksphaStandard = -dfc_data["ksphaStandard"];   # coefficients of dynamic measurment (up to 3rd order) of the standard method
ksphaNominal  = -dfc_data["ksphaNominal"];    # coefficients of nominal trajectory (1st order)

startStitched = dfc_data["startStitched"];
startStandard = dfc_data["startStandard"];
startNominal  = dfc_data["startNominal"];

tauStitched   = dfc_data["tauStitched"];      # synchronization delay
tauStandard   = dfc_data["tauStandard"];
tauNominal    = dfc_data["tauNominal"];

dt_adc        = dfc_data["dt_adc"];
datatime      = vec(dfc_data["datatime"]);

ksphaNominal  = InterpTrajTime(ksphaNominal , dt_Nominal , startNominal  + tauNominal  * dt_adc, datatime);
ksphaStandard = InterpTrajTime(ksphaStandard, dt_Standard, startStandard + tauStandard * dt_adc, datatime);
ksphaStitched = InterpTrajTime(ksphaStitched, dt_Stitched, startStitched + tauStitched * dt_adc, datatime);

bfieldNominal  = traj2grad(ksphaNominal , dt_adc; dim=1) / γ * 1e3; # convert to mT/m
bfieldStandard = traj2grad(ksphaStandard, dt_adc; dim=1) / γ * 1e3; 
bfieldStitched = traj2grad(ksphaStitched, dt_adc; dim=1) / γ * 1e3;

# 1. read the *.seq file
seq = read_seq(seq_file)[3:end-3]; # read_seq function from KomaMRI.jl, return a struct of Sequence
HighOrderMRI.KomaMRI.plot_seq(seq) # plot_seq function from KomaMRI.jl, plot the Sequence

nSample = size(datatime, 1);    
nTerm       = 16;

dfcNominal  = bfieldNominal'  * 1e-3; # mT => T
dfcStandard = bfieldStandard' * 1e-3; 
dfcStitched = bfieldStitched' * 1e-3; 

t = dt_adc * (nSample-1);

GR_dfcNominal  = reshape([KomaMRIBase.Grad(dfcNominal[idx,:] , t, dt_adc/2, dt_adc/2, 0) for idx=1:9], :, 1);
GR_dfcStandard = reshape([KomaMRIBase.Grad(dfcStandard[idx,:], t, dt_adc/2, dt_adc/2, 0) for idx=1:nTerm], :, 1);
GR_dfcStitched = reshape([KomaMRIBase.Grad(dfcStitched[idx,:], t, dt_adc/2, dt_adc/2, 0) for idx=1:nTerm], :, 1);

hoseqNominal  = HO_Sequence(seq);   # hoseq
hoseqStandard = HO_Sequence(seq);   
hoseqStitched = HO_Sequence(seq);   

# copy the 1st-order nominal gradient data from the seq object to the hoseq object
hoseqNominal.GR_dfc[2:4, :]  = hoseqNominal.SEQ.GR;     
hoseqStandard.GR_dfc[2:4, :] = hoseqStandard.SEQ.GR;    
hoseqStitched.GR_dfc[2:4, :] = hoseqStitched.SEQ.GR;    

# "5" is the index of the readout block in the spiral sequence
hoseqNominal.GR_dfc[1:9,5]  = GR_dfcNominal;   
hoseqStandard.GR_dfc[:,5] = GR_dfcStandard;             
hoseqStitched.GR_dfc[:,5] = GR_dfcStitched;            

plt_seq(hoseqNominal)
plt_seq(hoseqStandard)
plt_seq(hoseqStitched)

# kspha used for recon, extracted from the hoseq object
_, _, _, kspha_nominal  = get_kspace(hoseqNominal;  Δt=1);
_, _, _, kspha_standard = get_kspace(hoseqStandard; Δt=1);
_, _, _, kspha_stitched = get_kspace(hoseqStitched; Δt=1);


#########################################################################################
# 2. Simulation with dynamic field changes and ΔB0 map (quadratic-shaped)
#     a. setting the coil-sensitivity used in simulation
#     b. we use dynamic fields from field monitoring (stitching)
#     c. recon the signal with NuFFT immeadiately after simulation
#########################################################################################
T = Float64;
nX = nY = 200; nZ = 1;  # matrix size for recon
Δx = Δy = 1e-3; Δz = 2e-3;
# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 1        # set phantom down-sample factor to 5
location = 0.65
BHO = BlochHighOrder("1111", true, true)                       # ["1111"] turn on all order terms of dynamic field change. turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D", x=0.2, y=0.2, z=10.0) # setting for Phantom: decide which phantom file to use, loading phantom from src/phantom/mat folder                                
# setting the coil sensitivity used in the simulation
csm_type  = :real_32cha;      # a simulated birdcage coil-sensitivity
csm_nCoil = 32;              # 9-channel
csm_nRow  = 4;
csm_nCol  = 8;

db0_type  = :quadratic;     
db0_max   = :100.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map

# 1. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # set Δw to 0 if B0=false
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation
# plt_phantom(obj, :ρ)

# 2. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params();
sim_params["sim_method"]  = BHO;      # using "BlochHighOrder" for simulation with high-order terms
sim_params["return_type"] = "mat";    # setting with "mat", return the signal data for all channel
sim_params["precision"]   = "f64";
sim_params["gpu"]         = true;     # using GPU for simulation
sim_params["gpu_device"]  = 0;        # set the GPU device number, if using GPU for simulation
sim_params["Nblocks"]     = 2000;     # the number of blocks for GPU simulation, set according to the GPU memory

# 3. simulate
signal    = simulate(obj, hoseqStitched, sys; sim_params);          
raw       = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos   = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha   = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)

# 4. Adding noise to signal data
snr       = 10;
data      = AddNoise(signal[:,:,1], snr);
nSample, nCha = size(data);

# show the effect of noise
raw       = signal_to_raw_data(reshape(data, (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos   = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha   = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)

#########################################################################################
# 3. Reconstruction with HighOrderOp
#     a. obtain the ΔB₀ map
#     b. construct the encoding operator (HighOrderOp), defined in our HighOrderMRI.jl
#     c. reconstruct the signal with the encoding operator
#########################################################################################
# Coil-Sensitivity Map
coil    = csm_Real_32cha(217, 181, verbose=true);
csm     = get_center_crop(coil, nX, nY)[end:-1:1, :, :];
fig_csm = plt_images(abs.(csm); dim=3, nRow=csm_nRow, nCol=csm_nCol)

# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
b0map     = brain_phantom2D_reference(phantom, :Δw, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); 
                    location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
b0map     = rotl90(b0map);
fig_b0map = plt_B0map(-b0map)

# Proton-density map (reference)
x_ref   = brain_phantom2D_reference(phantom, :ρ, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); location=location, ss=ss);
x_ref   = rotl90(x_ref);
fig_ref = plt_image(x_ref)

#############################################################################
# HighOrderOp, the extended signal model for high-order terms
#############################################################################
nSample, nCha = size(data);
# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
datatime = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);
# Add a negative sign, as the phase term used in the model is positive.
kspha_stitched = -kspha_stitched;  
kspha_standard = -kspha_standard;
kspha_nominal  = -kspha_nominal;
b0             = -b0map;            

# For HighOrderOp
x, y = 1:nX, 1:nY;
x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- nX/2 .- 1, y .- nY/2 .- 1
x, y = x * Δx, y * Δy; 
gridding = Grid(nX=nX, nY=nY, nZ=nZ, Δx=Δx, Δy=Δy, Δz=Δz, x=T.(x), y=T.(y), z=T.(z));

recon_terms = "1111";  
#= The string "111" is a three-digit flag that indicates whether the 0th, 1st, and 2nd order terms of 
a measurement are used. For example, "110" means only the 0th and 1st order terms are used. =#
nBlock  = 20;   # the number is set according to the GPU memory.
use_gpu = true;
verbose = false;

# solver = "admm"; regularization = "TV"; iter = 5; λ = 1e-1;
solver = "cgnr"; regularization = "L2"; iter = 20; λ = 1e-9;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

weight = SampleDensity(kspha_stitched'[2:3,:], (nX, nY));

HOOp = HighOrderOp(gridding, T.(kspha_stitched'), T.(datatime); recon_terms=recon_terms, #k_nominal=kspha_nominal', 
                        nBlock=nBlock, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);


# recon with stitched measurement, with density weighting, with ΔB₀
@time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(x); vmaxp=99.9, title="Stitched")

# ssim, compared with the ground truth
@info "SSIM" SSIM=HO_SSIM(x_ref, abs.(x))

# to see the reconstructed images at each iteration
iter = 5;
plt_image(reshape(abs.(recParams[:solverInfo].x_iter[iter+1]), nX, nY); vmaxp=99.9, title="iter $(iter)")


# recon with Nominal trajectory, with density weighting, with ΔB₀
recon_terms = "0100";
HOOp = HighOrderOp(gridding, T.(kspha_nominal'), T.(datatime); recon_terms=recon_terms, 
                        nBlock=nBlock, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
@time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(x); vmaxp=99.9, title="Nominal")
