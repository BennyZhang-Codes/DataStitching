using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=false, T2=false, ss=11)                       # turn on B0, turn off T2, set phantom subsampling to 5
BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
csmtype= :real_32cha
nCoil   = 32; nrows=4; ncols=8;
maxOffresonance = 0.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map

fov = 150
grad_scale = 1 / 7.5
Nx = Ny = Int(round(500 * grad_scale, digits=0));
Sx = Sy = fov / Nx;

dir = "Figures/debug"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
hoseq_stitched = demo_hoseq(dfc_method=:Stitched, r=30)[4:end]   # :Stitched
hoseq_standard = demo_hoseq(dfc_method=:Standard, r=30)[4:end]   # :Standard
hoseq_stitched.SEQ.GR[1:2,5] = hoseq_stitched.SEQ.GR[1:2,5] * grad_scale
plot_seq(hoseq_stitched)

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
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
imgs_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
p_images = plot_imgs_subplots(abs.(imgs_nufft), nrows, ncols; title="$(nCoil) coils: NUFFT recon", height=400, width=800)
fig = plt_image(rotl90(sqrt.(sum(imgs_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54,title="")
fig.tight_layout(pad=0, w_pad=0, h_pad=0)
# fig.savefig("$(dir)/R$(R)_NUFFT.png", pad_inches=0, dpi=300, bbox_inches="tight")
# CSM 

coil = csm_Real_32cha(99, 83);
coil = get_center_crop(coil, Nx, Ny);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
p_smap_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")
# plt_images(permutedims(abs.(sensitivity[:,:,1,:]), (3, 1,2)))



x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(Sx,Sy));
plt_image(rotl90(x_ref ), title="reference ρ")


acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

T = Float64;
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = "L2"  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = 1e-4
recParams[:iterations] = 30
recParams[:solver] = "cgnr"
recParams[:toeplitz] = false
recParams[:oversamplingFactor] = 2
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
plt_image(rotl90(rec); title="HighOrderOp")
