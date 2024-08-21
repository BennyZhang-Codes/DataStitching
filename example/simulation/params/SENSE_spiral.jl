using KomaHighOrder
using MRIReco
using MAT
using PyPlot
import RegularizedLeastSquares: SolverInfo
using ImageDistances, ImageQualityIndexes
############################################################################################## 
# Setup
############################################################################################## 
R=4
simtype = SimType(B0=false, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
csmtype= :real_32cha
nCoil   = 32; nrows=4; ncols=8;
BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 0.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;
shape = (Nx, Ny);
T = Float64;
dir = "Figures/params_sense_0820/"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
seq = load_seq(seqname="demo", r=R)[4:end]
hoseq = HO_Sequence(seq)
plot_seq(hoseq)

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance, csmtype=csmtype, nCoil=nCoil)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"
sim_params["Nblocks"]    = 10000
# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
imgs_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(sqrt.(sum(imgs_nufft.^2; dims=3))[:,:,1]); title="Sim: $(BHO.name), R=$(R), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")


# CSM 
coil = csmtype == :real_32cha ? csm_Real_32cha(217, 181) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
coil = get_center_crop(coil, Nx, Ny);
sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
p_smap_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")

x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
f_ref = plt_image(rotl90(x_ref), title="Reference Image (Simulation)")

acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;



solver = "cgnr"
regularization = "L2"
λ = 1.e-2



#########################################################################
# Reconstruction with multiCoil and sense maps
#########################################################################
@info "Solver: $(solver), Regularization: $(regularization), λ: $(λ)"
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = 50
recParams[:solver] = solver
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true);

recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
plt_image(rotl90(rec); title="multiCoil Sense")


#########################################################################
# Reconstruction with multiCoil and sense maps, HighOrderOp
#########################################################################


