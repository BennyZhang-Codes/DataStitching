using HighOrderMRI
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances
import HighOrderMRI: get_grads
using Interpolations
using CUDA
gpu_idx = 0
CUDA.device!(gpu_idx)

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
ss = 1        # set phantom down-sample factor to 5
location = 0.8
BHO = BlochHighOrder("010", true, true)                          # ["111"] turn on all order terms of dynamic field change. turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.1, y=0.1, z=0.2) # setting for Phantom: decide which phantom file to use, loading phantom from src/phantom/mat folder                                
# setting the coil sensitivity used in the simulation
csm_type  = :fan      # coil type
csm_nCoil = 1         # number of partitions (split in fan shape)
csm_nRow  = nothing  
csm_nCol  = nothing  

db0_type  = :quadratic;     
db0_max   = :0.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map

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
sim_params["gpu_device"]  = gpu_idx;    

# 4. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);          
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
# fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)

# 5. Adding noise to signal data
snr = Inf;

kdata = AddNoise(signal[:,:,1], snr);
nSample, nCha = size(kdata);
# show the effect of noise
raw = signal_to_raw_data(reshape(kdata, (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
# fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)


#########################################################################################
# Prepare for reconstruction
#########################################################################################
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

csm = Array{Complex{T}}(undef, nX, nY, 1);
csm[:,:,1] = Complex{T}.(headmask);
csm[:,:,1] = Complex{T}.(ones(nX, nY));
#########################################################################################
# Reconstruction with HighOrderOp
#########################################################################################
nSample, nCha = size(kdata);
# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
# get nominal trajectory
seq = hoseqStitched.SEQ;
dt = hoseqStitched.SEQ.DEF["adc_DwellTime"];
k  = nothing
for s in seq
    if HighOrderMRI.is_ADC_on(s)
        t = collect(0:dt:s.DUR[1]);
        gx, gy, gz = get_grads(s, t)
        kx = cumtrapz(ones(length(t)-1)'*dt, gx')
        ky = cumtrapz(ones(length(t)-1)'*dt, gy')
        kz = cumtrapz(ones(length(t)-1)'*dt, gz')
        k = [kx' ky' kz'] * γ
    end
end

ksphaNominal = zeros(Float64, size(k, 1), 9);
ksphaNominal[:,2:4] = k[:,1:3];
startNominal = 0.;  # plt_kspha(ksphaNominal, dt)

_, _, _, ksphaStitched = get_kspace(hoseqStitched; Δt=1);
_, _, _, ksphaStandard = get_kspace(hoseqStandard; Δt=1);
datatime = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);
# Add a negative sign, as the phase term used in the model is positive.
ksphaStitched = -ksphaStitched;  
ksphaStandard = -ksphaStandard;
ksphaNominal  = -ksphaNominal;
b0            = -b0map;
# 
# plt_kspha(ksphaNominal, dt)
# plt_kspha(ksphaStitched, dt)
# plt_kspha_com(ksphaNominal[1:88000, :], ksphaStitched, dt)

# For HighOrderOp
x, y = 1:nX, 1:nY;
x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- nX/2 .- 1, y .- nY/2 .- 1
x, y = x * Δx, y * Δy; 
gridding = Grid(nX=nX, nY=nY, nZ=nZ, Δx=Δx, Δy=Δy, Δz=Δz, x=T.(x), y=T.(y), z=T.(z));

#= The string "111" is a three-digit flag that indicates whether the 0th, 1st, and 2nd order terms of 
a measurement are used. For example, "110" means only the 0th and 1st order terms are used. =#
Nblocks = 20;   # the number is set according to the GPU memory.
use_gpu = true;
verbose = true;


intermode = AkimaMonotonicInterpolation()
iter_max  = 20;
JumpFact  = 6;
Δτ_min    = T.(0.001);
λ         = T.(0);
BHO       = BlochHighOrder("010")

tauNominal  = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaNominal[:,1:9]),  T.(startNominal/dt),  T.(dt/dt);
            intermode=intermode, JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks, use_gpu=use_gpu, verbose=verbose)



csm = Array{Complex{T}}(undef, nX, nY, 1);
# csm[:,:,1] = Complex{T}.(headmask);
csm[:,:,1] = Complex{T}.(ones(nX, nY));
solver = "cgnr"; regularization = "L2"; iter = 20; λ = 0e-7;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)
weight   = SampleDensity(ksphaStitched'[2:3,:], (nX, nY));

HOOp = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BHO, 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);

# # recon with stitched measurement, with density weighting, with ΔB₀
@time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
plt_image(abs.(x); vmaxp=99.9, title="w/  ΔB₀, Stitched: 010, w/  density weighting")

HOOp = HighOrderOp(gridding, T.(ksphaNominal[1:88000, 1:9]'), T.(datatime); sim_method=BHO, 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);

# # recon with stitched measurement, with density weighting, with ΔB₀
@time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
plt_image(abs.(x); vmaxp=99.9, title="w/  ΔB₀, Nominal: 010, w/  density weighting")





results = Dict("tauNominal"=>tauNominal);
MAT.matwrite("$(outpath)/ksphaNominal.mat", results)
