using KomaHighOrder, MRISampling, MRIReco, MRICoilSensitivities
import KomaHighOrder.MRIBase: rawdata
R = 2
simtype  = SimType(B0=false, T2=false, ss=5)
coil_type= :real
Nparts   = 32; nrows=4; ncols=8;
Ncoils = Nparts
BHO_name = "000"
folder   = "spiral_$(coil_type)_Ncoils$(Ncoils)"
path     = "$(@__DIR__)/src/demo/demo_sense/results/$folder"
if ispath(path) == false mkdir(path) end

seq = demo_seq(seq="spiral", r=R)
# seq = demo_seq()
hoseq = HO_Sequence(seq)
Nx = 150
Ny = 150

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

coil_images = Array{Float32, 3}(undef, Nx, Ny, Ncoils);
signal = zeros(ComplexF64, sum(seq.ADC.N), Ncoils);
for coil_idx = 1:Ncoils
    @info "Simulating coil $(coil_idx)..."
    obj = brain_phantom2D(brain2D(); ss=simtype.ss, location=0.8, coil_type=coil_type, coil_idx=coil_idx, Nparts=Ncoils, overlap=0); 
    obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

    # simulate
    signal[:, coil_idx] = simulate(obj, hoseq, sys; sim_params);
    protocolName = "spiral_R$(R)_$(BHO_name)_nominal_$(Ncoils)-$(coil_idx)"
    raw = signal_to_raw_data(signal[:, coil_idx:coil_idx], hoseq, :nominal)
    coil_images[:,:,coil_idx] = reconstruct_2d_image(raw, Nx=Nx, Ny=Ny);

#     p = plot_image(coil_images[:,:,coil_idx]; title="$(protocolName)", height=400, width=450)
#     savefig(p,  "$(path)/$(protocolName).svg",format="svg", height=400, width=450)
end


raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "spiral_R$(R)_Ncoils$(Ncoils)"
raw.params["protocolName"] = filename
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)
p_coil_images = plot_imgs_subplots(abs.(coil_images), nrows, ncols; title="$(Ncoils) coils: NUFFT recon", height=400, width=800)
p_sos = plot_image(sqrt.(sum(coil_images.^2; dims=3))[:,:,1]; title="$(Ncoils) coils: NUFFT recon, SOS", height=400, width=450)
savefig(p_sos        , "$(path)/$(filename)-nufft_SOS.svg"        , format="svg", height=400, width=450)
savefig(p_coil_images, "$(path)/$(filename)-nufft.svg", format="svg", height=400, width=800)



raw = RawAcquisitionData(ISMRMRDFile("$(path)/spiral_R2_Ncoils$(Ncoils).mrd"));

r = 2
acqData = AcquisitionData(raw); # raw = RawAcquisitionData(mrd);
# acqData = convertUndersampledData(sample_kspace(acqData, r, "regular"))
acqData.traj[1].circular = false;
p_traj = plot_traj2d(acqData.traj[1]; height=400, width=400)
savefig(p_traj, "$(path)/$(filename)-traj.svg", format="svg", height=400, width=400)

shape = (Nx, Ny);
T = Float32;
#############################################################################
# recon with the coil sensitivities as the same used in the simulation
#############################################################################
coil = coil_type == :real ? coil = RealCoilSensitivity_32cha(217, 181) : coil = BirdcageSensitivity(217, 181, Ncoils, relative_radius=1.5);
coil = get_center_crop(coil, Nx, Ny);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, Ncoils);
for c = 1:Ncoils
    sensitivity[:,:,1,c] = coil[:,:,c]'
end

p_smap_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(Ncoils) coils: Coil Sensitivity (Simulation), SOS")
p_smap = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(Ncoils) coils: Coil Sensitivity (Simulation)")

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
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, Ncoils));
img_recon = Array{ComplexF32,2}(undef, Nx, Ny);
img_recon = reconstruction(acqData, params).data;
p_img_recon = plot_image(abs.(img_recon[:,:]), title="$(Ncoils) coils: Sense recon, R = $(r), sim-sensitivity")

savefig(p_smap     , "$(path)/$(filename)_simSmap.svg"            , format="svg", height=400, width=800)
savefig(p_smap_sos , "$(path)/$(filename)_simSmap_sos.svg"        , format="svg", height=400, width=450)
savefig(p_img_recon, "$(path)/$(filename)-simSmap_Sense.svg", format="svg", height=400, width=450)


#############################################################################
# recon with the coil sensitivities estimated from the raw data by espirit
#############################################################################
# espirit
raw_R1 = RawAcquisitionData(ISMRMRDFile("$(path)/spiral_R1_Ncoils$(Ncoils).mrd"));
# NumberOfSamples  = Int64(raw_R1.profiles[1].head.number_of_samples);
# NumberOfProfiles = Int64(length(raw_R1.profiles));
# NumberOfChannels = Ncoils;
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
p_smap_espirit_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(Ncoils) coils: Coil Sensitivity (espirit), SOS")
p_smap_espirit = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(Ncoils) coils: Coil Sensitivity (espirit)")

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
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, Ncoils));
img_recon = Array{ComplexF32,2}(undef, Nx, Ny);
img_recon = reconstruction(acqData, params).data;
p_img_recon_espirit = plot_image(abs.(img_recon[:,:]), title="$(Ncoils) coils: Sense recon, R = $(r), espirit-sensitivity")


savefig(p_smap_espirit     , "$(path)/$(filename)_espiritSmap.svg"            , format="svg", height=400, width=800)
savefig(p_smap_espirit_sos , "$(path)/$(filename)_espiritSmap_sos.svg"        , format="svg", height=400, width=450)
savefig(p_img_recon_espirit, "$(path)/$(filename)_espiritSmap_Sense.svg", format="svg", height=400, width=450)







seq = demo_seq(seq="spiral", r=1)
hoseq = HO_Sequence(seq)
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder("000");
# sim_params["sim_method"] = Bloch();
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

obj = brain_phantom2D(brain2D(); ss=3, location=0.5, coil_type=coil_type, coil_idx=1, Nparts=Ncoils); 
# obj = brain_phantom2D(brain2D(); ss=3);
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;  # TODO: fix the bug: gre 

# simulate
sig = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(sig, hoseq, :nominal)
# raw = signal_to_raw_data(sig, seq)
plot_signal(raw)
img = reconstruct_2d_image(raw; Nx=Nx, Ny=Ny);
p = plot_image(img)

# fft & plot kspace
kspace = KomaMRI.fftc(img); plot_image(abs.(kspace).^0.1)

# undersampling
sub = collect(range(1; stop=256, step=2))
plot_image(abs.(KomaMRI.ifftc(kspace[:, sub])))
plot_image(abs.(KomaMRI.fftc(KomaMRI.ifftc(kspace[:, sub]))).^0.1;)

# undersampling
acqData = AcquisitionData(raw)
acqData = convertUndersampledData(sample_kspace(acqData, 2, "regular"))
acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
# Nx, Ny = raw.params["reconSize"][1:2]
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
rec = reconstruction(acqData, recParams)

rec = reconstruction(a0, recParams)
image3d  = reshape(rec.data, Nx, Ny, :)
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
p = plot_image(image2d)
plot_image(abs.(KomaMRI.fftc(image2d)).^0.1;)




