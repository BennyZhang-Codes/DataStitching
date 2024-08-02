using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 200.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;

dir = "Figures/out"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
hoseq_stitched = demo_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = demo_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
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
img_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(img_nufft); title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
fig_b0map = plt_image(rotl90(B0map), title="B0map [-$maxOffresonance, $maxOffresonance] Hz")
x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));

acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_dfc_adc_stitched = get_kspace(hoseq_stitched; Δt=1);
_, _, _, K_dfc_adc_standard = get_kspace(hoseq_standard; Δt=1);

times = KomaMRIBase.get_adc_sampling_times(hoseq_stitched.SEQ);

tr_nominal          = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_stitched     = Trajectory(K_dfc_adc_stitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_standard     = Trajectory(K_dfc_adc_standard'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L1"  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = 1e-2
recParams[:iterations] = 90
recParams[:solver] = "fista"
recParams[:solverInfo] = SolverInfo(vec(ComplexF32.(x_ref)), store_solutions=true)

imgs = Array{Float32,3}(undef, 5, Nx, Ny);
# 1. ideally, with BlochHighOrder("000"), we use nominal trajectory and a null B0map
Op1 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map).*0, grid=1);
# 2. include ΔB₀, with BlochHighOrder("000"), we use nominal trajectory
Op2 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);
# 3. include both ΔB₀ and 0th/1st order terms of stitched DFC, with BlochHighOrder("110")
Op3 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);
# 4. include both ΔB₀ and all order terms of stitched DFC, with BlochHighOrder("111")
Op4 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);
# 5. include both ΔB₀ and all order terms of standard DFC, with BlochHighOrder("111")
Op5 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);
Ops = [Op1, Op2, Op3, Op4, Op5];
titles = ["HighOrderOp, stitched: 000, without ΔB₀",
          "HighOrderOp, stitched: 000, with ΔB₀",
          "HighOrderOp, stitched: 110 with ΔB₀",
          "HighOrderOp, stitched: 111 with ΔB₀",
          "HighOrderOp, standard: 111 with ΔB₀"]
for idx in eachindex(Ops)
    recParams[:encodingOps] = reshape([Ops[idx]], 1,1);
    @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
    imgs[idx, :, :] = rotl90(rec);
    plt_image(rotl90(rec); title=titles[idx])
end

MAT.matwrite(dir*"/recon_fully_sampled_fista_90_L1_1e-2.mat", Dict("imgs"=>imgs))
