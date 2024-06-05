using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISampling

simtype  = SimType(B0=false, T2=false, ss=5)
mask_type = :fan
overlap  = 0
Nparts   = 9
Npartsx  = 3
Npartsy  = 4
Ncoils = mask_type == :rect ? Npartsx * Npartsy : Nparts
BHO_name = "000"
hoseq    = demo_hoseq()
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_Ncoils$(Ncoils)"
if mask_type == :fan
    folder   = "$(mask_type)_Ncoils$(Ncoils)_overlap$(overlap)"
elseif mask_type == :rect
    folder   = "$(mask_type)_$(Npartsx)_$(Npartsy)"
end

path     = "$(@__DIR__)/src/demo/demo_sense/$folder"

raw_file = "$(path)/$(filename).mrd"
@assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
raw = RawAcquisitionData(ISMRMRDFile(raw_file));
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

# espirit
# acqDataCart = regrid(acqData, (Nx, Ny); cgnr_iter=3);
# sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98);

# fan mask
if mask_type == :fan
    mask = get_fan_mask(217, 181, Ncoils; overlap=overlap);
elseif mask_type == :rect
    mask = get_rect_mask(217, 181, Npartsx, Npartsy)
end
c1 = KomaHighOrder.get_center_range(217, Nx)
c2 = KomaHighOrder.get_center_range(181, Ny)
mask = mask[c1, c2, :]

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, Ncoils);
for c = 1:Ncoils
    sensitivity[:,:,1,c] = mask[:,:,c]'
end

p_mask = plot_imgs_subplots(mask, 3,3)

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

p_img_cg = plot_image(abs.(img_cg[:,:]))

savefig(p_mask,  "$(path)/$(raw.params["protocolName"])-Mask.svg", format="svg", height=400, width=450)
savefig(p_img_cg,  "$(path)/$(raw.params["protocolName"])-Recon.svg", format="svg", height=400, width=450)




