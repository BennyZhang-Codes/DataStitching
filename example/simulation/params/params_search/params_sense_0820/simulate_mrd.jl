using KomaHighOrder
using MAT
using PyPlot
##############################################################################################
# Setup
##############################################################################################
R = 4
simtype = SimType(B0=false, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
csmtype = :real_32cha
nCoil   = 32; nrows=4; ncols=8;
BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 200.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;

dir = "Figures/out"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
seq = load_seq(seqname="demo", r=R)
hoseq = HO_Sequence(seq)[4:end]
plot_seq(hoseq)
# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, csmtype=csmtype, nCoil=nCoil, B0type=:quadratic, maxOffresonance=maxOffresonance)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);


fig_nufft = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]); title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")
