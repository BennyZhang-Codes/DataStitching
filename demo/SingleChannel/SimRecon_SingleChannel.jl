using KomaHighOrder, MRIReco, MAT
import RegularizedLeastSquares: SolverInfo

#########################################################################################
# 1. Load sequence and dynamic field data
#     a. Load sequence file (*.seq format from pulseq)
#     b. Load dynamic field data (*.mat format)
#     c. Create a HO_Sequence object as defined in KomaHighOrder.jl by combining the 
#        sequence and dynamic field data.
#########################################################################################

seq_file = "$(@__DIR__)/demo/SingleChannel/1mm_R1.seq"   # *.seq file is the pulseq's sequence file
# under-sampling factor: R = 1
# inplane resolution: 1 mm x 1 mm
# FOV: 150 mm x 150 mm
# readout duration: 88 ms

dfc_file = "$(@__DIR__)/demo/SingleChannel/1mm_R1.mat"   # *.mat file contains the dynamic field data from both stitching method and the standard method.
# The dynamic field data is stored in the *.mat file with the following keys:
#= 
"dt"               [s], time interval between two time points.
"gradStandard"     [mT], first order term: x, y, from the standard method, complex-valued.
"gradStitched"     [mT], first order term: x, y, from the stitched method, complex-valued.
"gradSegStitched"  [mT], first order term: x, y, from the stitched method, complex-valued. A matrix, one column for each segment.
"ntStandard"       number of time points of the standard method.
"ntStitched"       number of time points of the stitched method. A vector, one element for each segment.
"skopeStandard"    [mT], fields of dynamic measurment (up to 2nd order) of the standard method.
"skopeStitched"    [mT], fields of dynamic measurment (up to 2nd order) of the stitched method.
=#

print("$(@__DIR__)")
# 1. read the *.seq file and reverse the sign of the gradient (axis x)
seq = read_seq(seq_file)[4:end];   # read_seq function from KomaMRI.jl, return a struct of Sequence
seq.GR[1,:] = -seq.GR[1,:];        # reverse the sign of the gradient (axis x)
plot_seq(seq)                      # plot_seq function from KomaMRI.jl, plot the Sequence


grad = MAT.matread(dfc_file);
Δt = grad["dt"];
nGradSample, nTerm = size(grad["skopeStitched"])              
dfcStitched = grad["skopeStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
dfcStandard = grad["skopeStandard"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
ntStitched  = grad["ntStitched"]
ntStandard  = grad["ntStandard"]
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
#     a. we use dynamic fields from field monitoring (stitching)
#     b. recon the signal with NuFFT immeadiately after simulation
#########################################################################################
# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
BHO = BlochHighOrder("111", true, true)                          # ["111"] turn on all order terms of dynamic field change. turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # setting for Phantom: decide which phantom file to use, loading phantom from src/phantom/mat folder
maxOffresonance = 50.                                            # set the maximum off-resonance frequency in Hz for quadratic B0 map

# 1. sequence
plot_seq(hoseqStitched)
Nx=Ny=150;     # matrix size for recon

# 2. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params();
sim_params["sim_method"]  = BHO;      # using "BlochHighOrder" for simulation with high-order terms
sim_params["return_type"] = "mat";    # setting with "mat", return the signal data for all channel
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);      
fig_nufft = plt_image(rotl90(img_nufft))


#########################################################################################
# 3. Reconstruction with HighOrderOp
#     a. obtain the ΔB₀ map
#     b. construct the encoding operator (HighOrderOp), defined in our KomaHighOrder.jl
#     c. reconstruct the signal with the encoding operator
#########################################################################################
# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
B0map = brain_phantom2D_reference(phantom; ss=ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
fig_b0map = plt_image(rotl90(B0map))
# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom; ss=ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));
fig_ref = plt_image(rotl90(x_ref))

acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
_, K_nominal_adc, _, K_dfc_adc_stitched = get_kspace(hoseqStitched; Δt=1);
times = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);

# create the trajectories for the nominal and the stitching measurement
tr_nominal      = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_stitched = Trajectory(K_dfc_adc_stitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

# Construct the encoding operator (HighOrderOp)
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);

# Reconstruction parameters
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = 1.e-3
recParams[:iterations] = 20
recParams[:solver] = "cgnr"  # "cgnr", "admm"
recParams[:solverInfo] = SolverInfo(vec(ComplexF32.(x_ref)), store_solutions=true);
recParams[:encodingOps] = reshape([Op], 1,1);

# Run reconstruction
@time img_HOOP = abs.(reconstruction(acqData, recParams).data[:,:]);
plt_image(rotl90(img_HOOP); title="w/  ΔB₀, stitched: 111")