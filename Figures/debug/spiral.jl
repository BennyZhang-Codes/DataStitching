using KomaHighOrder, MRISampling, MRIReco, MRICoilSensitivities
import KomaHighOrder.MRIBase: rawdata
R = 2
simtype  = SimType(B0=false, T2=false, ss=5)
csmtype= :real_32cha
nCoil   = 32; nrows=4; ncols=8;
BHO_name = "000"
folder   = "spiral_$(csmtype)_nCoil$(nCoil)"
path     = "$(@__DIR__)/demo/demo_sense/results/$folder"
if ispath(path) == false mkpath(path) end

seq = load_seq(seqname="spiral", r=R)
hoseq = HO_Sequence(seq)

hoseq = demo_hoseq(dfc_method=:Stitched, r=30)[4:end]   # :Stitched
hoseq.SEQ.GR[1:2,5] = hoseq.SEQ.GR[1:2,5] * 1/7.5
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
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)
images = recon_2d(raw, Nx=Nx, Ny=Ny);

p_images = plot_imgs_subplots(abs.(images), nrows, ncols; title="$(nCoil) coils: NUFFT recon", height=400, width=800)
p_sos = plot_image(sqrt.(sum(images.^2; dims=3))[:,:,1]; title="$(nCoil) coils: NUFFT recon, SOS", height=400, width=450)
savefig(p_sos   , "$(path)/$(filename)-nufft_SOS.svg", format="svg", height=400, width=450)
savefig(p_images, "$(path)/$(filename)-nufft.svg"    , format="svg", height=400, width=800)



raw = RawAcquisitionData(ISMRMRDFile("$(path)/spiral_R2_nCoil$(nCoil).mrd"));

r = 2
acqData = AcquisitionData(raw); # raw = RawAcquisitionData(mrd);
acqData.traj[1].circular = false;
p_traj = plot_traj2d(acqData.traj[1]; height=400, width=400)
savefig(p_traj, "$(path)/$(filename)-traj.svg", format="svg", height=400, width=400)


Nx, Ny, _ = raw.params["reconSize"]
shape = (Nx, Ny);
T = Float32;
#############################################################################
# recon with the coil sensitivities as the same used in the simulation
#############################################################################
coil = csmtype == :real_32cha ? csm_Real_32cha(217, 181) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
coil = get_center_crop(coil, 150, 150);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = imresize(transpose(coil[:,:,c]), shape)
end

p_smap_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(nCoil) coils: Coil Sensitivity (Simulation), SOS")
p_smap_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")
p_smap_pha = plot_imgs_subplots(angle.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 50
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 1
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
img_recon = Array{ComplexF32,2}(undef, Nx, Ny);
img_recon = reconstruction(acqData, params).data;
p_img_recon = plot_image(abs.(img_recon[:,:]), title="$(nCoil) coils: Sense recon, R = $(r), sim-sensitivity")

savefig(p_smap_mag , "$(path)/$(filename)_simSmap_mag.svg"  , format="svg", height=400, width=800)
savefig(p_smap_pha , "$(path)/$(filename)_simSmap_pha.svg"  , format="svg", height=400, width=800)
savefig(p_smap_sos , "$(path)/$(filename)_simSmap_sos.svg"  , format="svg", height=400, width=450)
savefig(p_img_recon, "$(path)/$(filename)-simSmap_Sense.svg", format="svg", height=400, width=450)


#############################################################################
# recon with the coil sensitivities estimated from the raw data by espirit
#############################################################################
# espirit
raw_R1 = RawAcquisitionData(ISMRMRDFile("$(path)/spiral_R1_nCoil$(nCoil).mrd"));
# NumberOfSamples  = Int64(raw_R1.profiles[1].head.number_of_samples);
# NumberOfProfiles = Int64(length(raw_R1.profiles));
# NumberOfChannels = nCoil;
# raw_R1.params["trajectory"] = "custom";
# kdata = rawdata(raw_R1);
# _, ktraj_adc = get_kspace(demo_seq(seq="spiral", r=1));
# tr           = Trajectory(T.(ktraj_adc[:,1:2]'), NumberOfProfiles, NumberOfSamples, circular=false, cartesian=false);
# dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
# dat[1,1,1]   = reshape(kdata,:,NumberOfChannels);
# a            = AcquisitionData(tr, dat, encodingSize=(Nx,Ny));
a = AcquisitionData(raw_R1)
a = regrid(a, (Nx,Ny); cgnr_iter=3);
sensitivity = espirit(a, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.95);
p_smap_espirit_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(nCoil) coils: Coil Sensitivity (espirit), SOS")
p_smap_espirit_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (espirit)")
p_smap_espirit_pha = plot_imgs_subplots(angle.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (espirit)")

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 50
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 1
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));
img_recon = Array{ComplexF32,2}(undef, Nx, Ny);
img_recon = reconstruction(acqData, params).data;
p_img_recon_espirit = plot_image(abs.(img_recon[:,:]), title="$(nCoil) coils: Sense recon, R = $(r), espirit-sensitivity")


savefig(p_smap_espirit_mag , "$(path)/$(filename)_espiritSmap_mag.svg"  , format="svg", height=400, width=800)
savefig(p_smap_espirit_pha , "$(path)/$(filename)_espiritSmap_pha.svg"  , format="svg", height=400, width=800)
savefig(p_smap_espirit_sos , "$(path)/$(filename)_espiritSmap_sos.svg"  , format="svg", height=400, width=450)
savefig(p_img_recon_espirit, "$(path)/$(filename)_espiritSmap_Sense.svg", format="svg", height=400, width=450)







