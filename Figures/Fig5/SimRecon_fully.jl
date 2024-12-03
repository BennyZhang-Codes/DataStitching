using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances

outpath = "$(@__DIR__)/Figures/Fig5/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory
############################################################################################## 
# Setup
############################################################################################## 
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
location = 0.8;
Nx = Ny = 150;

# setting the coil sensitivity used in the simulation
csm_type  = :fan;      # a simulated birdcage coil-sensitivity
csm_nCoil = 1;         # 1-channel
csm_nRow  = 1;
csm_nCol  = 1;

db0_type  = :quadratic;     
db0_max   = :200.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map


solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;

# 1. sequence
hoseq_stitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["gpu"] = false;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseq_stitched, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(img_nufft); title="Sim: $(BHO.name), Δw: [-$db0_max,$db0_max] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(db0_max)_reconNUFFT.svg", width=550,height=500,format="svg")

# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (1., 1.); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
fig_b0map = plt_B0map(rotl90(B0map))

# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=ss);
fig_ref = plt_image(rotl90(x_ref))

# 5. reconstruction
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
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
# recParams[:solverInfo] = SolverInfo(vec(ComplexF32.(x_ref)), store_solutions=true);


Op = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=20, fieldmap=Matrix(B0map).*0, grid=1, use_gpu=false, verbose=true);
recParams[:encodingOps] = reshape([Op], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
plt_image(rotl90(rec))



Nblocks=9;
##### w/o ΔB₀
# 1. nominal trajectory, BlochHighOrder("000")
Op1 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 2. stitched trajectory, BlochHighOrder("110")
Op2 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 3. stitched trajectory, BlochHighOrder("111")
Op3 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);
# 4. standard trajectory, BlochHighOrder("111")
Op4 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map).*0, grid=1);

##### with ΔB₀
# 5. nominal trajectory, BlochHighOrder("000")
Op5 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 6. stitched trajectory, BlochHighOrder("110")
Op6 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 7. stitched trajectory, BlochHighOrder("111")
Op7 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
# 8. standard trajectory, BlochHighOrder("111")
Op8 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1);
Ops = [Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8];

imgs = Array{Float32,3}(undef, length(Ops), Nx, Ny);
labels = ["w/o ΔB₀, stitched: 000",
          "w/o ΔB₀, stitched: 110",
          "w/o ΔB₀, stitched: 111",
          "w/o ΔB₀, standard: 111",
          "w/  ΔB₀, stitched: 000",
          "w/  ΔB₀, stitched: 110",
          "w/  ΔB₀, stitched: 111",
          "w/  ΔB₀, standard: 111",];
for idx in eachindex(Ops)
    recParams[:encodingOps] = reshape([Ops[idx]], 1,1);
    @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
    imgs[idx, :, :] = rotl90(rec);
    plt_image(rotl90(rec); title=labels[idx])
end

MAT.matwrite("$(outpath)/fully_$(solver)_$(iter)_$(regularization)_$(λ).mat", Dict("imgs"=>imgs, "labels"=>labels))
