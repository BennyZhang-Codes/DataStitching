using KomaHighOrder
using MRIReco, MRICoilSensitivities, MRISampling

Ncoils = 30
BHO_name = "000"
hoseq = demo_hoseq()
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_Ncoils$(Ncoils)"

folder = "Ncoils$Ncoils"
path = "$(@__DIR__)/src/demo/demo_sense/$folder"

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
mask = get_fan_mask(Nx, Ny, Ncoils)';
p_mask = plot_image(mask)
sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, Ncoils);
for c = 1:Ncoils
    sensitivity[:,:,1,c] = ((mask .== c) .* 1)
end



# s = [sensitivity[:,:,1,i] for i=1:Ncoils];
# s = [s[1] s[2] s[3]];
# plot_image(abs.(rotr90(s)); title="sensitivity - espirit")
# plot_image(abs.(((mask .== 1) .* 1)'))

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




