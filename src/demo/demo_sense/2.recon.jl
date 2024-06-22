using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISampling

simtype  = SimType(B0=false, T2=false, ss=5)
coil_type= :birdcage
overlap  = 0
Nparts   = 9; nrows=3; ncols=3;
Npartsx  = 3
Npartsy  = 4
Ncoils = coil_type == :rect ? Npartsx * Npartsy : Nparts
BHO_name = "000"
hoseq    = demo_hoseq()
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_Ncoils$(Ncoils)"
if coil_type == :fan
    folder   = "$(coil_type)_Ncoils$(Ncoils)_overlap$(overlap)"
elseif coil_type == :rect
    folder   = "$(coil_type)_$(Npartsx)_$(Npartsy)"
elseif coil_type == :birdcage
    folder   = "$(coil_type)_Ncoils$(Ncoils)"
end

path     = "$(@__DIR__)/src/demo/demo_sense/$folder"

raw_file = "$(path)/$(filename).mrd"
@assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
raw = RawAcquisitionData(ISMRMRDFile(raw_file));
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);



# fan mask
if coil_type == :fan
    coil = get_fan_mask(217, 181, Ncoils; overlap=overlap);
elseif coil_type == :rect
    coil = get_rect_mask(217, 181, Npartsx, Npartsy)
elseif coil_type == :birdcage
    coil = BirdcageSensitivity(217, 181, Ncoils)
end
c1 = KomaHighOrder.get_center_range(217, Nx)
c2 = KomaHighOrder.get_center_range(181, Ny)
coil = coil[c1, c2, :]

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, Ncoils);
for c = 1:Ncoils
    sensitivity[:,:,1,c] = coil[:,:,c]'
end

p_smap = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(Ncoils) coils: Coil Sensitivity", height=400, width=450)

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
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, Ncoils));

img_cg = Array{ComplexF32,2}(undef, Nx, Ny);

img_cg = reconstruction(acqData, params).data;

p_img_cg = plot_image(abs.(img_cg[:,:]), title="$(Ncoils) coils: Sense-type Recon", height=400, width=450)

savefig(p_smap,  "$(path)/$(raw.params["protocolName"])-CoilSens.svg", format="svg", height=400, width=450)
savefig(p_img_cg,  "$(path)/$(raw.params["protocolName"])-Recon.svg", format="svg", height=400, width=450)





# espirit
acqDataCart = regrid(acqData, (Nx, Ny); cgnr_iter=3);
sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.97);

p_smap_espirit = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(Ncoils) coils: Coil Sensitivity (espirit)", height=400, width=450)

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
params[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, Ncoils));

img_cg = Array{ComplexF32,2}(undef, Nx, Ny);

img_cg = reconstruction(acqData, params).data;

p_img_cg_espirit = plot_image(abs.(img_cg[:,:]), title="$(Ncoils) coils: Sense-type Recon (espirit)", height=400, width=450)

savefig(p_smap_espirit,  "$(path)/$(raw.params["protocolName"])-CoilSens_espirit.svg", format="svg", height=400, width=450)
savefig(p_img_cg_espirit,  "$(path)/$(raw.params["protocolName"])-Recon_espirit.svg", format="svg", height=400, width=450)




