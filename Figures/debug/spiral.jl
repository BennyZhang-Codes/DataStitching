using KomaHighOrder, MRIReco
import KomaHighOrder.MRIBase: rawdata
R = 30
simtype  = SimType(B0=false, T2=false, ss=5)
csmtype= :real_32cha
nCoil   = 32; nrows=4; ncols=8;
BHO_name = "000"
folder   = "spiral_$(csmtype)_nCoil$(nCoil)"
path     = "Figures/debug"; if ispath(path) == false mkpath(path) end     # output directory
if ispath(path) == false mkpath(path) end

seq = load_seq(seqname="demo", r=R)
hoseq = HO_Sequence(seq)

# hoseq = demo_hoseq(dfc_method=:Stitched, r=R)[4:end]   # :Stitched
hoseq.SEQ.GR[1:2,8] = hoseq.SEQ.GR[1:2,8] * 1/7.5
plot_seq(hoseq)

Nx = 150
Ny = 150

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

obj = brain_hophantom2D(BrainPhantom(); ss=simtype.ss, location=0.8, csmtype=csmtype, nCoil=nCoil); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

# simulate
signal = simulate(obj, hoseq, sys; sim_params);

raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "spiral_R$(R)_nCoil$(nCoil)"
raw.params["protocolName"] = filename
# mrd = ISMRMRDFile("$(path)/$(filename).mrd")
# save(mrd, raw)
images = recon_2d(raw, Nx=Nx, Ny=Ny);

p_images = plot_imgs_subplots(abs.(images), nrows, ncols; title="$(nCoil) coils: NUFFT recon", height=400, width=800)
p_sos = plot_image(sqrt.(sum(images.^2; dims=3))[:,:,1]; title="$(nCoil) coils: NUFFT recon, SOS", height=400, width=450)
# PlotlyJS.savefig(p_sos   , "$(path)/$(filename)-nufft_SOS.svg", format="svg", height=400, width=450)
# PlotlyJS.savefig(p_images, "$(path)/$(filename)-nufft.svg"    , format="svg", height=400, width=800)



# raw = RawAcquisitionData(ISMRMRDFile("$(path)/spiral_R2_nCoil$(nCoil).mrd"));

# r = 2
acqData = AcquisitionData(raw); # raw = RawAcquisitionData(mrd);
acqData.traj[1].circular = false;
p_traj = plot_traj2d(acqData.traj[1]; height=400, width=400)
# savefig(p_traj, "$(path)/$(filename)-traj.svg", format="svg", height=400, width=400)


# Nx, Ny, _ = raw.params["reconSize"]
shape = (Nx, Ny);
T = Float32;
#############################################################################
# recon with the coil sensitivities as the same used in the simulation
#############################################################################
coil = csmtype == :real_32cha ? csm_Real_32cha(217, 181) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
coil = get_center_crop(coil, Nx, Ny);

# coil = csmtype == :real_32cha ? csm_Real_32cha(1085, 905) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
# coil = get_center_crop(coil, 466, 480);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end

p_smap_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-4)
params[:iterations] =30
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 2
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
img_recon = Array{ComplexF32,2}(undef, Nx, Ny);
img_recon = reconstruction(acqData, params).data;
p_img_recon = plot_image(abs.(img_recon[:,:]), title="$(nCoil) coils: Sense recon, R = $(R), sim-sensitivity")


