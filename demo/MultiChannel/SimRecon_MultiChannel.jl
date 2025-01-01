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
nX = nY = 200; nZ = 1;  # matrix size for recon
Δx = Δy = 1e-3; Δz = 2e-3;
# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
location = 0.65
BHO = BlochHighOrder("111", true, true)                          # ["111"] turn on all order terms of dynamic field change. turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D", x=0.2, y=0.2, z=0.2) # setting for Phantom: decide which phantom file to use, loading phantom from src/phantom/mat folder                                
# setting the coil sensitivity used in the simulation
csm_type  = :birdcage;      # a simulated birdcage coil-sensitivity
csm_nCoil = 16;              # 9-channel
csm_nRow  = 4;
csm_nCol  = 4;

db0_type  = :quadratic;     
db0_max   = :100.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map

# 1. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # set Δw to 0 if B0=false
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation
# plot_phantom_map(obj, :ρ)
# plot_phantom_map_csm(obj, :mag; coil_idx=1)  # [:pha, :real, :imag] plot the magnitude of coil-sensitivity in channel 1.

# 2. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params();
sim_params["sim_method"]  = BHO;      # using "BlochHighOrder" for simulation with high-order terms
sim_params["return_type"] = "mat";    # setting with "mat", return the signal data for all channel
sim_params["precision"]   = "f64";
sim_params["Nblocks"] = 1000;
sim_params["gpu_device"]  = 0;    

# 3. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);          
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)

# 4. Adding noise to signal data
snr = 15;
data = signal[:,:,1];
nSample, nCha = size(data);

signalAmpl = sum(abs.(data), dims=1)/ nSample;
data = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));

# show the effect of noise
raw = signal_to_raw_data(reshape(data, (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)
#########################################################################################
# 3. Reconstruction with HighOrderOp
#     a. obtain the ΔB₀ map
#     b. construct the encoding operator (HighOrderOp), defined in our KomaHighOrder.jl
#     c. reconstruct the signal with the encoding operator
#########################################################################################
# Coil-Sensitivity Map
coil = csm_Birdcage(217, 181, csm_nCoil, verbose=true);
csm  = get_center_crop(coil, nX, nY)[end:-1:1, :, :];
fig_csm = plt_images(abs.(csm); dim=3, nRow=csm_nRow, nCol=csm_nCol)

# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
b0map = brain_phantom2D_reference(phantom, :Δw, (200., 200.), (1., 1.); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
b0map = rotl90(b0map);
fig_b0map = plt_B0map(b0map)

# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom, :ρ, (200., 200.), (1., 1.); location=location, ss=ss);
fig_ref = plt_image(rotl90(x_ref))

#############################################################################
# HighOrderOp, the extended signal model for high-order terms
#############################################################################
# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
_, kspha_nominal, _, kspha_stitched = get_kspace(hoseqStitched; Δt=1);
_,             _, _, kspha_standard = get_kspace(hoseqStandard; Δt=1);
datatime = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);
# Add a negative sign, as the phase term used in the model is positive.
kspha_stitched = -kspha_stitched;  
kspha_standard = -kspha_standard;
kspha_nominal  = -kspha_nominal;
b0             = -b0map;            

nSample, nCha = size(data);

# For HighOrderOp
T = Float64;

x, y = 1:nX, 1:nY;
x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- nX/2 .- 1, y .- nY/2 .- 1


x, y = 1:nX, 1:nY;
x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- nX/2 .- 1, y .- nY/2 .- 1
x, y = x * Δx, y * Δy; 
gridding = Grid(nX=nX, nY=nY, nZ=nZ, Δx=Δx, Δy=Δy, Δz=Δz, x=T.(x), y=T.(y), z=T.(z));

BHO = BlochHighOrder("111");  
#= The string "111" is a three-digit flag that indicates whether the 0th, 1st, and 2nd order terms of 
a measurement are used. For example, "110" means only the 0th and 1st order terms are used. =#
Nblocks = 20;   # the number is set according to the GPU memory.
use_gpu = true;
verbose = true;

solver = "admm"; regularization = "TV"; iter = 5; λ = 1e-15;
solver = "cgnr"; regularization = "L2"; iter = 5; λ = 0*1e-9;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
# recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=false)

weight = SampleDensity(kspha_stitched'[2:3,:], (nX, nY));

HOOp = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);

# recon with stitched measurement, with density weighting, with ΔB₀
@time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(x1); vmaxp=99.9, title="w/  ΔB₀, stitched: 111, w/  density weighting")

solverinfo = recParams[:solverInfo];
solverinfo.

# recon with stitched measurement, without density weighting, with ΔB₀
@time x2 = recon_HOOp(HOOp, Complex{T}.(data), recParams);
plt_image(abs.(x2); vmaxp=99.9, title="w/  ΔB₀, stitched: 111, w/o  density weighting")


# recon with stitched measurement, with density weighting, without ΔB₀
HOOp = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0.*0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
@time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(x); vmaxp=99.9, title="w/o ΔB₀, stitched: 111, w/  density weighting")


# recon with nominal trajectory, with density weighting, with ΔB₀
BHO = BlochHighOrder("000");  # "000" indicates that no measured field dynamics are used, only nominal kspace trajectory is used.
HOOp = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
@time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(x); vmaxp=99.9, title="w/  ΔB₀, nominal, w/  density weighting")
# you can try different recons like: "011", "101", "110"...





Nblocks=10;
##### with ΔB₀
# 1. nominal trajectory, BlochHighOrder("000")
Op1 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 2. stitched trajectory, BlochHighOrder("110")
Op2 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 3. stitched trajectory, BlochHighOrder("111")
Op3 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 4. standard trajectory, BlochHighOrder("111")
Op4 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStandard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);


##### w/o ΔB₀
# 5. nominal trajectory, BlochHighOrder("000")
Op5 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 6. stitched trajectory, BlochHighOrder("110")
Op6 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 7. stitched trajectory, BlochHighOrder("111")
Op7 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 8. standard trajectory, BlochHighOrder("111")
Op8 = HighOrderOp((nX, nY), tr_nominal, tr_ksphaStandard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);

Ops = [Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8];


imgs = Array{T,3}(undef, length(Ops), nX, nY);
labels = [ "wB0_nominal",  "wB0_stitched_110",  "wB0_stitched_111",  "wB0_standard_111",
          "woB0_nominal", "woB0_stitched_110", "woB0_stitched_111", "woB0_standard_111",];
for idx in eachindex(Ops)
    Op = DiagOp(Ops[idx], numChan) ∘ S 
    recParams[:encodingOps] = reshape([Op], 1,1);
    @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
    imgs[idx, :, :] = rotl90(rec);
    plt_image(rotl90(rec); title=labels[idx])
end



#############################################################################
# Save the results
#############################################################################
# ΔB₀ map
B0map = brain_phantom2D_reference(phantom, :Δw, (200., 200.), (1., 1.); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
B0map = rotl90(B0map);

x_ref = brain_phantom2D_reference(phantom, :ρ, (200., 200.), (1., 1.); location=location, ss=ss);
x_ref = rotl90(x_ref);

headmask = brain_phantom2D_reference(phantom, :headmask, (200., 200.), (1., 1.); location=location, ss=ss);
headmask = rotl90(headmask);

MAT.matwrite("$(@__DIR__)/demo/MultiChannel/snr$(snr)_$(solver)_$(iter)_$(regularization)_$(λ).mat", 
            Dict("imgs"=>imgs, "labels"=>labels, "csm"=>csm, "signal"=>data, "B0map"=>B0map, "x_ref"=>x_ref, "headmask"=>headmask))
