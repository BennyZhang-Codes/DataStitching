using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances
using CUDA
CUDA.device!(1)

outpath = "$(@__DIR__)/Figures/Fig5/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory
############################################################################################## 
# Setup
############################################################################################## 
T = Float64;
nX = nY = 150; nZ = 1;  # matrix size for recon
Δx = Δy = 1e-3; Δz = 2e-3;
# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
location = 0.8
BHO = BlochHighOrder("111", true, true)                          # ["111"] turn on all order terms of dynamic field change. turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # setting for Phantom: decide which phantom file to use, loading phantom from src/phantom/mat folder                                
# setting the coil sensitivity used in the simulation
csm_type  = :birdcage;      # a simulated birdcage coil-sensitivity
csm_nCoil = 16;              # 9-channel
csm_nRow  = 4;
csm_nCol  = 4;

db0_type  = :quadratic;     
db0_max   = :150.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map


solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;

# 1. sequence
hoseqStitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseqStandard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # set Δw to 0 if B0=false
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params();
sim_params["sim_method"]  = BHO;      # using "BlochHighOrder" for simulation with high-order terms
sim_params["return_type"] = "mat";    # setting with "mat", return the signal data for all channel
sim_params["precision"]   = "f64";
sim_params["Nblocks"] = 1000;
sim_params["gpu_device"]  = 1;    

# 4. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);          
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)

# 5. Adding noise to signal data
snr = 10;

kdata = AddNoise(signal[:,:,1], snr);
nSample, nCha = size(kdata);
# show the effect of noise
raw = signal_to_raw_data(reshape(kdata, (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)


#########################################################################################
# Prepare for reconstruction
#########################################################################################
# Coil-Sensitivity Map
coil = csm_Birdcage(217, 181, csm_nCoil, verbose=true);
csm  = get_center_crop(coil, nX, nY)[end:-1:1, :, :];
fig_csm = plt_images(abs.(csm); dim=3, nRow=csm_nRow, nCol=csm_nCol)

# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
b0map = brain_phantom2D_reference(phantom, :Δw, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
b0map = rotl90(b0map);
fig_b0map = plt_B0map(b0map)

# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom, :ρ, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); location=location, ss=ss);
x_ref = rotl90(x_ref);
fig_ref = plt_image(x_ref)

headmask = brain_phantom2D_reference(phantom, :headmask, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); location=location, ss=ss);
headmask = rotl90(headmask);
fig_headmask = plt_image(headmask)

#########################################################################################
# Reconstruction with HighOrderOp
#########################################################################################
nSample, nCha = size(kdata);
# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
_, kspha_nominal, _, kspha_stitched = get_kspace(hoseqStitched; Δt=1);
_,             _, _, kspha_standard = get_kspace(hoseqStandard; Δt=1);
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
weight   = SampleDensity(kspha_stitched'[2:3,:], (nX, nY));
BHO      = BlochHighOrder("111");  
#= The string "111" is a three-digit flag that indicates whether the 0th, 1st, and 2nd order terms of 
a measurement are used. For example, "110" means only the 0th and 1st order terms are used. =#
Nblocks = 20;   # the number is set according to the GPU memory.
use_gpu = true;
verbose = true;

solver = "admm"; regularization = "TV"; iter = 20; λ = 5e-1;
# solver = "cgnr"; regularization = "L2"; iter = 20; λ = 1e-7;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

HOOp = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);

# # recon with stitched measurement, with density weighting, with ΔB₀
# @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
# plt_image(abs.(x); vmaxp=99.9, title="w/  ΔB₀, stitched: 111, w/  density weighting")

# recon with stitched measurement, with density weighting, with ΔB₀
@time x = recon_HOOp(HOOp, Complex{T}.(kdata), recParams);
plt_image(abs.(x); vmaxp=99.9, title="w/  ΔB₀, stitched: 111, w/  density weighting")

# ssim, compared with the ground truth
println("SSIM: ", HO_SSIM(x_ref, abs.(x)))



##### with ΔB₀
# 1. nominal trajectory, BlochHighOrder("000")
HOOp1 = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("000"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# 2. stitched trajectory, BlochHighOrder("110")
HOOp2 = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("110"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# 3. stitched trajectory, BlochHighOrder("111")
HOOp3 = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("111"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# 4. standard trajectory, BlochHighOrder("111")
HOOp4 = HighOrderOp(gridding, T.(kspha_standard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("111"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
##### w/o ΔB₀
# 5. nominal trajectory, BlochHighOrder("000")
HOOp5 = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("000"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0.*0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# 6. stitched trajectory, BlochHighOrder("110")
HOOp6 = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("110"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0.*0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# 7. stitched trajectory, BlochHighOrder("111")
HOOp7 = HighOrderOp(gridding, T.(kspha_stitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("111"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0.*0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# 8. standard trajectory, BlochHighOrder("111")
HOOp8 = HighOrderOp(gridding, T.(kspha_standard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("111"), tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0.*0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);

Ops = [HOOp1, HOOp2, HOOp3, HOOp4, HOOp5, HOOp6, HOOp7, HOOp8];


imgs = Array{Complex{T},3}(undef, nX, nY, length(Ops));
labels = [ "wB0_nominal",  "wB0_stitched_110",  "wB0_stitched_111",  "wB0_standard_111",
          "woB0_nominal", "woB0_stitched_110", "woB0_stitched_111", "woB0_standard_111",];

solver = "cgnr"; regularization = "L2"; iter = 20; λ = 1e-7;
solver = "admm"; regularization = "TV"; iter = 20; λ = 5e-1;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

for idx in eachindex(Ops)
    @info "Recon with $(labels[idx])"
    HOOp = Ops[idx]
    # @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
    @time x = recon_HOOp(HOOp, Complex{T}.(kdata), recParams);
    imgs[:, :, idx] = x;
    plt_image(abs.(x); vmaxp=99.9, title=labels[idx])
end

results = Dict("imgs"=>imgs, "labels"=>labels, "csm"=>csm, "b0"=>b0, "datatime"=>datatime, "snr"=>snr, 
            "solver"=>solver, "iter"=>iter, "regularization"=>regularization, "lambda"=>λ,
            "ksphaNominal"=>kspha_nominal, "ksphaStandard"=>kspha_standard, "ksphaStitched"=>kspha_stitched);

MAT.matwrite("$(outpath)/1p0_150_r1_snr$(snr)_$(solver)_$(iter)_$(regularization)_$(λ).mat", results)
