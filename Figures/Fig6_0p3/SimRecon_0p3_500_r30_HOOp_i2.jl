using HighOrderMRI
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances

outpath = "$(@__DIR__)/Figures/Fig6/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory
############################################################################################## 
# Setup
############################################################################################## 
T = Float64;
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 3        # set phantom down-sample factor to 3
location = 0.8;                                              
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.1, y=0.1, z=0.2) # decide which phantom file to use
Nx = Ny = 500;
# setting the coil sensitivity used in the simulation
csm_type  = :gaussian_grid_block;      
# csm_nCoil = 256;             
# csm_nRow  = 16;
# csm_nCol  = 16;
# csm_nBlock = 4;
# csm_radius = 5;

csm_nCoil = 100;             
csm_nRow  = 10;
csm_nCol  = 10;
csm_nBlock = 4;
csm_radius = 4.5;
csm_gpu   = true;


db0_type  = :quadratic;     
db0_max   = :100.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map


#############################################################################
# 1. Simulation
#############################################################################
# 1. sequence
hoseq_stitched = load_hoseq(seqname="spiral", r=30, dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = load_hoseq(seqname="spiral", r=30, dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
@time obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, csm_radius=csm_radius, csm_nBlock=csm_nBlock, csm_gpu=csm_gpu,
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"]   = "f64"
sim_params["Nblocks"] = 200;
sim_params["Nthreads"] = 1;

# 4. simulate
signal = simulate(obj, hoseq_stitched, sys; sim_params);
data = signal[:,:,1];
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
# fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
# fig_cha = plt_images(mapslices(rotl90, img_nufft, dims=[1,2]); dim=3, nRow=csm_nRow, nCol=csm_nCol)


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
coil = csm_Gaussian_grid_block(724, 604, csm_nCoil, true; nRow=csm_nRow, nCol=csm_nCol, nBlock=csm_nBlock, relative_radius=csm_radius, verbose=true);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = reshape(permutedims(coil, (2,1,3)), Nx, Ny, 1, csm_nCoil);
csm = mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]);
fig_csm = plt_images(csm; dim=3, nRow=csm_nRow, nCol=csm_nCol)

smap = permutedims(sensitivity, [1,2,4,3])[:,:,:,1];# (nY, nX, nCha, 1)
yik_sos = sum(abs.(conj(smap) .* img_nufft); dims=3)[:,:,1]; # coil combine
fig_coilcombine = plt_image(rotl90(abs.(yik_sos)))

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (0.3,0.3); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
plt_B0map(rotl90(B0map))
x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (0.3,0.3); location=location, ss=ss);
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
Δx = Δy = 0.3e-3;

solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:solverInfo] = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);
recParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));

# Nblocks = 20
# Op = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Δx=Δx, Δy=Δy, 
#                         Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map, grid=1, use_gpu=true, verbose=true);

# @time x = recon_HOOp(Op, acqData, recParams)
# plt_image(rotl90(abs.(x)))

# solverinfo = recParams[:solverInfo];
# plt_image(rotl90(reshape(abs.(solverinfo.x_iter[3]), Nx, Ny)))






Nblocks=10;
##### w/  ΔB₀
# 1. nominal trajectory, BlochHighOrder("000")
Op1 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map, grid=1, verbose=true);
# 2. stitched trajectory, BlochHighOrder("110")
Op2 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map, grid=1, verbose=true);
# 3. stitched trajectory, BlochHighOrder("111")
Op3 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map, grid=1, verbose=true);
# 4. standard trajectory, BlochHighOrder("111")
Op4 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map, grid=1, verbose=true);

##### w/o ΔB₀
# 5. nominal trajectory, BlochHighOrder("000")
Op5 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map.*0, grid=1, verbose=true);
# 6. stitched trajectory, BlochHighOrder("110")
Op6 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map.*0, grid=1, verbose=true);
# 7. stitched trajectory, BlochHighOrder("111")
Op7 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map.*0, grid=1, verbose=true);
# 8. standard trajectory, BlochHighOrder("111")
Op8 = HighOrderOp_i2((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(sensitivity[:,:,1,:]), fieldmap=B0map.*0, grid=1, verbose=true);
Ops = [Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8];


imgs = Array{T,3}(undef, length(Ops), Nx, Ny);
labels = [ "wB0_nominal",  "wB0_stitched_110",  "wB0_stitched_111",  "wB0_standard_111",
          "woB0_nominal", "woB0_stitched_110", "woB0_stitched_111", "woB0_standard_111",];

for idx in 1:8 #eachindex(Ops)
    @time rec = abs.(recon_HOOp(Ops[idx], acqData, recParams));
    imgs[idx, :, :] = rotl90(rec);
    plt_image(rotl90(rec); title=labels[idx])
end


HO_img_scale(rotl90(x_ref), imgs[3, :, :])
plt_image(abs.(imgs[3, :, :] - imgs[4, :, :])*10, vmax=maximum(imgs[3,:,:]),vmin=0)
plt_image(abs.(imgs[3, :, :] - rotl90(x_ref)), vmax=maximum(imgs[3,:,:]),vmin=0)

plt_image(abs.(HO_img_scale(rotl90(x_ref), imgs[4, :, :]) - rotl90(x_ref)))

HO_SSIM(imgs[4, :, :], rotl90(x_ref))

fig,ax = plt.subplots(1,1, figsize=(4,4))
ax.plot(tr_dfc_standard.nodes[2,:])
ax.plot(tr_dfc_stitched.nodes[2,:])
plt_image(abs.(imgs[3, :, :] - imgs[4, :, :]))

B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (0.3,0.3); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (0.3,0.3); location=location, ss=ss);
headmask = brain_phantom2D_reference(phantom, :headmask, (150., 150.), (0.3,0.3); location=location, ss=ss);
MAT.matwrite("$(outpath)/R30_snr$(snr)_$(solver)_$(iter)_$(regularization)_$(λ).mat", Dict("imgs"=>imgs, "labels"=>labels, "csm"=>csm, "signal"=>data, "B0map"=>rotl90(B0map), "x_ref"=>rotl90(x_ref), "headmask"=>rotl90(headmask)));



solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));

imgs_admm = Array{T,3}(undef, length(Ops), Nx, Ny);

idx = 3;
@time rec = abs.(recon_HOOp(Ops[idx], acqData, recParams));
imgs[idx, :, :] = rotl90(rec);
plt_image(rotl90(rec); title=labels[idx])

idx = 4;
@time rec = abs.(recon_HOOp(Ops[idx], acqData, recParams));
imgs[idx, :, :] = rotl90(rec);
plt_image(rotl90(rec); title=labels[idx])


# # ΔB₀ map
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (0.3,0.3); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
# B0map = rotl90(B0map);

x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (0.3,0.3); location=location, ss=ss);
# x_ref = rotl90(x_ref);

headmask = brain_phantom2D_reference(phantom, :headmask, (150., 150.), (0.3,0.3); location=location, ss=ss);
# headmask = rotl90(headmask);

MAT.matwrite("$(outpath)/R30_snr$(snr)_$(solver)_$(iter)_$(regularization)_$(λ).mat", Dict("imgs"=>imgs, "labels"=>labels, "csm"=>csm, "signal"=>data, "B0map"=>rotl90(B0map), "x_ref"=>rotl90(x_ref), "headmask"=>rotl90(headmask)));

