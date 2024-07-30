using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISimulation
using MAT, ImageQualityIndexes, ImageDistances
# using ProgressMeter

import ImageTransformations: imresize
############################################################################################## 
# Setup
############################################################################################## 
key = :Stitched
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use


dir = "Figures/out"; if ispath(dir) == false mkpath(dir) end
maxOffresonance = 200.
Nx = Ny = 150;


# 1. hoseq
hoseq = demo_hoseq(dfc_method=key)[4:end]   # :Standard or :Stitched

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation


# 4. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"
# 5. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img = recon_2d(raw);
p_image = plot_image(img; darkmode=true, title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")
