using KomaHighOrder, MRISampling, MRIReco

R = 1
simtype  = SimType(B0=false, T2=true, ss=3)
coil_type= :birdcage
Nparts   = 3; nrows=1; ncols=3;
Ncoils = Nparts
BHO_name = "000"
folder   = "gre_$(coil_type)_Ncoils$(Ncoils)"
path     = "$(@__DIR__)/src/demo/demo_sense/results/$folder"
if ispath(path) == false mkdir(path) end

seq = demo_seq(seq="gre", r=R)
# seq = demo_seq()
hoseq = HO_Sequence(seq)
Nx = seq.DEF["Nx"]
Ny = seq.DEF["Ny"] 

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

coil_images = Array{Float32, 3}(undef, Nx, Ny, Ncoils);
signal = zeros(ComplexF64, sum(seq.ADC.N), Ncoils);
for coil_idx = 1:Ncoils
    obj = brain_phantom2D(brain2D(); ss=simtype.ss, location=0.5, coil_type=coil_type, coil_idx=coil_idx, Nparts=Ncoils); 
    obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; # TODO: fix the bug: gre 

    # simulate
    signal[:, coil_idx] = simulate(obj, hoseq, sys; sim_params);
    protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_$(Ncoils)-$(coil_idx)"
    raw = signal_to_raw_data(signal[:, coil_idx:coil_idx], hoseq, :nominal)
    coil_images[:,:,coil_idx] = reconstruct_2d_image(raw, Nx=Nx, Ny=Ny);

    p = plot_image(coil_images[:,:,coil_idx]; title="$(protocolName)", height=400, width=450)
    savefig(p,  "$(path)/$(protocolName).svg",format="svg", height=400, width=450)
end


raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_Ncoils$(Ncoils)"
raw.params["protocolName"] = filename
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)
p_coil_images = plot_imgs_subplots(abs.(coil_images), nrows, ncols; title="$(Ncoils) coils: NUFFT recon", height=400, width=450)
p_sos = plot_img(sqrt.(sum(coil_images.^2; dims=3))[:,:,1]; title="SOS, $(name)", width=450, height=420)
savefig(p_sos, "$(path)/$(raw.params["protocolName"])-nufft_SOS.svg", format="svg", height=400, width=450)
savefig(p_coil_images,  "$(path)/$(raw.params["protocolName"])-nufft_multi-coils.svg", format="svg", height=400, width=450)



coil = BirdcageSensitivity(362, 302, Ncoils, relative_radius=1.5) # TODO: fix the bug: when encoding matrix is bigger than phantom matrix

c1 = KomaHighOrder.get_center_range(362, Nx)
c2 = KomaHighOrder.get_center_range(302, Ny)
coil = coil[c1, c2, :]

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, Ncoils);
for c = 1:Ncoils
    sensitivity[:,:,1,c] = coil[:,:,c]'
end
plot_img(abs.(sqrt.(sum(sensitivity.^2, dims=4)))[:,:,1,1])

p_smap = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(Ncoils) coils: Coil Sensitivity", height=400, width=450)


@info "reference reco"
acqData = AcquisitionData(raw);
# acqData = convertUndersampledData(sample_kspace(acqData, 2, "regular"))
acqData.traj[1].circular = false;
shape = (Nx, Ny);
T = Float32;
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 100
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 1
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, Ncoils));

img_cg = Array{ComplexF32,2}(undef, Nx, Ny);

img_cg = reconstruction(acqData, params).data;

p_img_cg = plot_image(abs.(img_cg[:,:]), title="$(Ncoils) coils: Sense-type Recon")

savefig(p_smap,  "$(path)/$(raw.params["protocolName"])-CoilSens.svg", format="svg", height=400, width=450)
savefig(p_img_cg,  "$(path)/$(raw.params["protocolName"])-Recon.svg", format="svg", height=400, width=450)





seq = demo_seq(seq="gre", r=1)
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





# different
a1 = AcquisitionData(raw)
a2 = convertUndersampledData(sample_kspace(a1, 2, "regular"))