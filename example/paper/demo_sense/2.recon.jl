using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISampling

simtype  = SimType(B0=false, T2=false, ss=5)
csmtype= :fan
overlap  = 0
nCoil   = 4; nrows, ncols = get_factors(nCoil);

BHO_name = "000"
if csmtype == :fan
    folder   = "$(csmtype)_nCoil$(nCoil)_overlap$(overlap)"
else
    folder   = "$(csmtype)_nCoil$(nCoil)"
end
path     = "$(@__DIR__)/demo/demo_sense/results/$folder"
if ispath(path) == false mkpath(path) end

filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_nCoil$(nCoil)"
raw_file = "$(path)/$(filename).mrd"
@assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
raw = RawAcquisitionData(ISMRMRDFile(raw_file));
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);



# csm
coil = load_csm(csmtype, 217, 181, nCoil; overlap=overlap, relative_radius=1.5)
coil = get_center_crop(coil, Nx, Ny);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
p_smap_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(nCoil) coils: Coil Sensitivity, SOS")
p_smap = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity")#, height=400, width=450)

@info "reference reco"
T = Float32;
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 1
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 2
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));

img_cg = Array{ComplexF32,2}(undef, Nx, Ny);
img_cg = reconstruction(acqData, params).data;
p_img_cg = plot_image(abs.(img_cg[:,:]), title="$(nCoil) coils: Sense-type Recon", height=400, width=450)

savefig(p_smap,  "$(path)/$(raw.params["protocolName"])-CoilSens.svg", format="svg", height=400, width=800)
savefig(p_img_cg,  "$(path)/$(raw.params["protocolName"])-Recon.svg", format="svg", height=400, width=450)
savefig(p_smap_sos,  "$(path)/$(raw.params["protocolName"])-CoilSens_sos.svg", format="svg", height=400, width=450)




# espirit
acqDataCart = regrid(acqData, (Nx, Ny); cgnr_iter=3);
sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.95);
p_smap_espirit_sos = plot_image(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]); title="$(nCoil) coils: Coil Sensitivity (espirit), SOS")
p_smap_espirit = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (espirit)")

@info "reference reco"
T = Float32;
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = (Nx, Ny)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 1
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 2
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCoil));

img_cg = Array{ComplexF32,2}(undef, Nx, Ny);

img_cg = reconstruction(acqData, params).data;

p_img_cg_espirit = plot_image(abs.(img_cg[:,:]), title="$(nCoil) coils: Sense-type Recon (espirit)", height=400, width=450)

savefig(p_smap_espirit,  "$(path)/$(raw.params["protocolName"])-CoilSens_espirit.svg", format="svg", height=400, width=800)
savefig(p_img_cg_espirit,  "$(path)/$(raw.params["protocolName"])-Recon_espirit.svg", format="svg", height=400, width=450)
savefig(p_smap_espirit_sos,  "$(path)/$(raw.params["protocolName"])-CoilSens_espirit_sos.svg", format="svg", height=400, width=450)




