using KomaHighOrder, MRIReco
import KomaHighOrder.MRIBase: rawdata
import RegularizedLeastSquares: SolverInfo
path     = "Figures/debug/out/R10_reduced_SENSE_CS_rectCSM"; if ispath(path) == false mkpath(path) end     # output directory

##############################################################################################
# Setup
##############################################################################################
T = Float32;
R = 30
matrix_origin = 500
matrix_target = 150
grad_scale    = 1/(matrix_origin/matrix_target)
Nx = Ny = matrix_target
shape = (Nx, Ny);

simtype  = SimType(B0=false, T2=false, ss=1)
csmtype= :rect; nCoil = 400; nrows=20; ncols=20;
maxOffresonance = 0.   

BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use

fileprefix = "spiral_R$(R)_nCoil$(nCoil)_R20"


##### 1. sequence
seq = load_seq(seqname="demo", r=R)
hoseq = HO_Sequence(seq)

# hoseq = demo_hoseq(dfc_method=:Stitched, r=R)[4:end]   # :Stitched
hoseq.SEQ.GR[1:2,8] = hoseq.SEQ.GR[1:2,8] * grad_scale
plot_seq(hoseq)

_, k_nominal, _, _ = get_kspace(hoseq; Δt=1);
fig_traj = plt_traj(k_nominal'; color_label="#CCCCCC")
fig_traj.savefig("$(path)/$(fileprefix)-traj.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)

##### 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, csmtype=csmtype, nCoil=nCoil, B0type=:quadratic, maxOffresonance=maxOffresonance); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

##### 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BHO;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["precision"]   = "f64"

##### 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);


raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);

fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54)
fig_cha = plt_images(permutedims(mapslices(rotl90, img_nufft,dims=[1,2]), [3, 1, 2]),width=8, height=8)

fig_sos.savefig("$(path)/$(fileprefix)-nufft_SOS.png", dpi=300, bbox_inches="tight", pad_inches=0)
fig_cha.savefig("$(path)/$(fileprefix)-nufft.png"    , dpi=300, bbox_inches="tight", pad_inches=0)



##### 5. plot trajectory
acqData = AcquisitionData(raw); # raw = RawAcquisitionData(mrd);
acqData.traj[1].circular = false;


#############################################################################
# recon with the coil sensitivities as the same used in the simulation
#############################################################################
# coil = csmtype == :real_32cha ? csm_Real_32cha(217, 181) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
coil = csm_Rect_binary(217, 181, nCoil, verbose=true);
coil = get_center_crop(coil, Nx, Ny);


sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end

fig_csm = plt_images(permutedims(mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]), [3, 1, 2]),width=8, height=8)
fig_csm.savefig("$(path)/$(fileprefix)-csm.png"    , dpi=300, bbox_inches="tight", pad_inches=0)

# x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(300, 300), target_resolution=(1,1));

solver = "admm"
reg = "TV"
iter = 50
λ = 1e-4
LSParams = Dict{Symbol,Any}()
LSParams[:reconSize]          = (Nx, Ny)
LSParams[:densityWeighting]   = true
LSParams[:reco]               = "multiCoil"
LSParams[:regularization]     = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
LSParams[:λ]                  = λ
LSParams[:iterations]         = iter
LSParams[:solver]             = solver
LSParams[:relTol]             = 0.0
LSParams[:oversamplingFactor] = 2
LSParams[:toeplitz]           = false
LSParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
# LSParams[:solverInfo]         = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);
@time rec = abs.(reconstruction(acqData, LSParams).data[:,:]);
fig = plt_image(rotl90(rec))

fig.savefig("$(path)/$(fileprefix)_SENSE_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)
# solverinfo = LSParams[:solverInfo];
# x_iters = solverinfo.x_iter[2:end];
# size(x_iters)


# L1-Wavelet regularized CS reconstruction
cs_solver = "admm"
cs_reg    = "L1"
cs_sparse = "Wavelet"
cs_λ      = 1.e-3
cs_iter   = 200

# cs_solver = "fista"
# cs_reg    = "L1"
# cs_sparse = "Wavelet"
# cs_λ      = 1.e-6
# cs_iter   = 100


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
CSParams[:senseMaps]          = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
CSParams[:normalizeReg]       = true
CSParams[:solverInfo]         = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);
@time rec = abs.(reconstruction(acqData, CSParams).data[:,:]);
fig = plt_image(rotl90(rec))


fig.savefig("$(path)/$(fileprefix)_CS_$(cs_solver)_$(cs_reg)_$(cs_sparse)_$(cs_iter)_$(cs_λ).png", dpi=300, bbox_inches="tight", pad_inches=0)
solverinfo = CSParams[:solverInfo];
x_iters = solverinfo.x_iter[2:end];
size(x_iters)

