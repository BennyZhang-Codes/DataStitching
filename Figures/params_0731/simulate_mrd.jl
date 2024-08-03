using KomaHighOrder
using MAT
using PyPlot
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

# hoseq_stitched.SEQ.GR[:,5] = hoseq_stitched.SEQ.GR[:,5] * 2
# hoseq_stitched.GR_dfc[:,5] = hoseq_stitched.GR_dfc[:,5] * 2
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
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);
fig_nufft = plt_image(rotl90(img_nufft); title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")
