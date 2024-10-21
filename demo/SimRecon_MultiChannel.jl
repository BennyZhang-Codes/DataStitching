using KomaHighOrder, MRIReco, MAT
import RegularizedLeastSquares: SolverInfo

#########################################################################################
# 1. Load sequence and dynamic field data
#     a. Load sequence file (*.seq format from pulseq)
#     b. Load dynamic field data (*.mat format)
#     c. Create a HO_Sequence object as defined in KomaHighOrder.jl by combining the 
#        sequence and dynamic field data.
#########################################################################################

seq_file = "$(@__DIR__)/demo/Sequence/1mm_R1.seq"        # *.seq file is the pulseq's sequence file
dfc_file = "$(@__DIR__)/demo/DynamicFields/1mm_R1.mat"   # *.mat file contains the dynamic field data from both stitching method and the standard method.
# The dynamic field data is stored in the *.mat file with the following keys:
#= 
"dt"               time interval between two time points [s].
"gradStandard"     first order term: x, y, from the standard method, complex-valued.
"gradStitched"     first order term: x, y, from the stitched method, complex-valued.
"gradSegStitched"  first order term: x, y, from the stitched method, complex-valued. A matrix, one column for each segment.
"ntStandard"       number of time points of the standard method.
"ntStitched"       number of time points of the stitched method. A vector, one element for each segment.
"skopeStandard"    dynamic measurment (up to 2nd order) of the standard method.
"skopeStitched"    dynamic measurment (up to 2nd order) of the stitched method.
=#

print("$(@__DIR__)")
# 1. read the *.seq file and reverse the sign of the gradient (axis x)
seq = read_seq(seq_file)[4:end];   # read_seq function from KomaMRI.jl, return a struct of Sequence
seq.GR[1,:] = -seq.GR[1,:];        # reverse the sign of the gradient (axis x)
plot_seq(seq)                      # plot_seq function from KomaMRI.jl, plot the Sequence


grad = MAT.matread(dfc_file);
Δt = grad["dt"];
nGradPoint, nTerm = size(grad["skopeStitched"])              
dfcStitched = grad["skopeStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
dfcStandard = grad["skopeStandard"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
ntStitched  = grad["ntStitched"]
ntStandard  = grad["ntStandard"]
t = Δt * (nGradPoint-1);
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
#     a. we use dynamic fields from stitching method for simulation
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
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);      
fig_nufft = plt_image(rotl90(img_nufft))