using KomaHighOrder, MRIReco
import KomaHighOrder.MRIBase: rawdata
import RegularizedLeastSquares: SolverInfo
outpath = "$(@__DIR__)/workplace/Parallel_Imaging/debug/out/R30_reduced_SENSE_gaussian_grid_pha_lin"; if ispath(outpath) == false mkpath(outpath) end     # output directory
##############################################################################################
# Setup
##############################################################################################
T = Float64;
R = 30
matrix_origin = 500
matrix_target = 500
grad_scale    = 1/(matrix_origin/matrix_target)
Nx = Ny = matrix_target
shape = (Nx, Ny);

simtype  = SimType(B0=false, T2=false, ss=3)
BHO      = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom  = BrainPhantom(prefix="brain3D724", x=0.1, y=0.1, z=0.2) # decide which phantom file to use
location = 0.8

# settings for phantom
csm_type  = :gaussian_grid;      # a simulated birdcage coil-sensitivity
csm_nCoil = 1600;              # 8-channel
csm_nRow  = 40;
csm_nCol  = 40;
csm_radius = 1.5;

db0_type  = :quadratic;     
db0_max   = :0.;

# brain_phantom2D_reference(phantom, :csm, (150., 150.), (1.,1.); location=location, ss=simtype.ss,
#                         csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
#                         db0_type=db0_type, db0_max=db0_max);


fileprefix = "spiral_radius$(csm_radius)_nCoil$(csm_nCoil)"

##### 1. sequence
seq = load_seq(seqname="spiral", r=R)
hoseq = HO_Sequence(seq)

# hoseq = demo_hoseq(dfc_method=:Stitched, r=R)[4:end]   # :Stitched
hoseq.SEQ.GR[1:2,8] = hoseq.SEQ.GR[1:2,8] * grad_scale
# plot_seq(hoseq)

_, k_nominal, _, _ = get_kspace(hoseq; Δt=1);
fig_traj = plt_traj(k_nominal'; color_label="#CCCCCC")
fig_traj.savefig("$(outpath)/Traj.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)

##### 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, csm_radius=csm_radius, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2πobj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

##### 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"]   = "f64"
sim_params["Nblocks"] = 1000;
##### 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);

fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
#fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow=csm_nRow, nCol=csm_nCol)

coil = csm_Gaussian_grid(724, 604, csm_nCoil; nRow=csm_nRow, nCol=csm_nCol, relative_radius=csm_radius, verbose=true);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = reshape(permutedims(coil, (2,1,3)), Nx, Ny, 1, csm_nCoil);
#fig_csm = plt_images(mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]); dim=3, nRow=csm_nRow, nCol=csm_nCol)

smap = permutedims(sensitivity, [1,2,4,3])[:,:,:,1];# (nY, nX, nCha, 1)
yik_sos = sum(abs.(conj(smap) .* img_nufft); dims=3)[:,:,1]; # coil combine
fig_coilcombine = plt_image(rotl90(abs.(yik_sos)))

#fig_coilcombine.savefig("$(outpath)/$(fileprefix)-nufft_coilcombine.png", dpi=300, bbox_inches="tight", pad_inches=0)
#fig_csm.savefig("$(outpath)/$(fileprefix)-csm.png", dpi=300, bbox_inches="tight", pad_inches=0)

#fig_sos.savefig("$(outpath)/$(fileprefix)-nufft_SOS.png", dpi=300, bbox_inches="tight", pad_inches=0)
#fig_cha.savefig("$(outpath)/$(fileprefix)-nufft.png"    , dpi=300, bbox_inches="tight", pad_inches=0)


#############################################################################
# recon with the coil sensitivities as the same used in the simulation
#############################################################################
acqData = AcquisitionData(raw, BHO; sim_params=sim_params); # raw = RawAcquisitionData(mrd);
acqData.traj[1].circular = false;

# x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (0.3,0.3); location=location, ss=simtype.ss);
# fig_ref = plt_image(rotl90(x_ref))
# fig_ref.savefig("$(outpath)/Reference.png", dpi=300, bbox_inches="tight", pad_inches=0)

# solver = "admm"; regularization = "TV"; iter = 20; λ = 1e-4
solver = "cgnr"; regularization = "L2"; iter = 20; λ = 1e-3
LSParams = Dict{Symbol,Any}()
LSParams[:reconSize]          = (Nx, Ny)
LSParams[:densityWeighting]   = true
LSParams[:reco]               = "multiCoil"
LSParams[:regularization]     = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
LSParams[:λ]                  = λ
LSParams[:iterations]         = iter
LSParams[:solver]             = solver
LSParams[:relTol]             = 0.0
LSParams[:oversamplingFactor] = 2
LSParams[:toeplitz]           = false
LSParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));
# LSParams[:solverInfo]         = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);
@time rec = abs.(reconstruction(acqData, LSParams).data[:,:]);
fig = plt_image(rotl90(rec))

fig.savefig("$(outpath)/$(fileprefix)_SENSE_$(solver)_$(regularization)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)
# solverinfo = LSParams[:solverInfo];
# x_iters = solverinfo.x_iter[2:end];
# size(x_iters)


# L1-Wavelet regularized CS reconstruction
cs_solver = "admm" ; cs_reg = "L1"; cs_sparse = "Wavelet"; cs_λ = 1.e-3; cs_iter = 200
cs_solver = "fista"; cs_reg = "L1"; cs_sparse = "Wavelet"; cs_λ = 1.e-3; cs_iter = 100


CSParams = Dict{Symbol, Any}()
CSParams[:reconSize]          = (Nx, Ny)
CSParams[:densityWeighting]   = true
CSParams[:reco]               = "multiCoil"
CSParams[:solver]             = cs_solver
CSParams[:regularization]     = cs_reg
CSParams[:sparseTrafo]        = cs_sparse
CSParams[:λ]                  = cs_λ
CSParams[:iterations]         = cs_iter
# CSParams[:ρ]                  = 0.1
CSParams[:absTol]             = 1.e-15
CSParams[:relTol]             = 1.e-9
CSParams[:tolInner]           = 1.e-9
CSParams[:adaptRho]           = true
CSParams[:iterationsInner]    = 100
CSParams[:oversamplingFactor] = 2
CSParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));
CSParams[:normalizeReg]       = true
CSParams[:solverInfo]         = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);
@time rec = abs.(reconstruction(acqData, CSParams).data[:,:]);
fig = plt_image(rotl90(rec))


fig.savefig("$(outpath)/$(fileprefix)_CS_$(cs_solver)_$(cs_reg)_$(cs_sparse)_$(cs_iter)_$(cs_λ).png", dpi=300, bbox_inches="tight", pad_inches=0)
solverinfo = CSParams[:solverInfo];
x_iters = solverinfo.x_iter[2:end];
size(x_iters)

