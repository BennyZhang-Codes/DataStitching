using KomaHighOrder

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

hoseq = demo_hoseq()

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["gpu"] = true;
sim_params["return_type"]="mat";
sim_params["Nblocks"] = 20;

Nx = Ny = 150
coil_images = Array{Float32, 3}(undef, Nx, Ny, nCoil);
signal = zeros(ComplexF64, sum(hoseq.SEQ.ADC.N), nCoil);
for coil_idx = 1:nCoil
    obj = brain_phantom2D(BrainPhantom(); ss=simtype.ss, location=0.8, csmtype=csmtype, coil_idx=coil_idx, nCoil=nCoil, overlap=overlap); 
    obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; 

    # simulate
    signal[:, coil_idx] = simulate(obj, hoseq, sys; sim_params);
    protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_$(nCoil)-$(coil_idx)"
    coil_images[:,:,coil_idx] = recon_2d(signal_to_raw_data(signal[:, coil_idx:coil_idx], hoseq, :nominal));

    p = plot_image(coil_images[:,:,coil_idx]; title="$(protocolName)", height=400, width=450)
    savefig(p,  "$(path)/$(protocolName).svg",format="svg", height=400, width=450)
end

raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_nCoil$(nCoil)"
raw.params["protocolName"] = filename
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)


p_sos = plot_image(abs.(sqrt.(sum(coil_images.^2; dims=3))[:,:,1]); title="$(nCoil) coils: NUFFT recon, SOS")
savefig(p_sos,  "$(path)/$(raw.params["protocolName"])-nufft_multi-coils_sos.svg", format="svg", height=400, width=450)
p_coil_images = plot_imgs_subplots(abs.(coil_images), nrows, ncols; title="$(nCoil) coils: NUFFT recon", height=400, width=450)
# savefig(p_coil_images,  "$(path)/$(raw.params["protocolName"])-nufft_multi-coils.svg", format="svg", height=400, width=450)
savefig(p_coil_images,  "$(path)/$(raw.params["protocolName"])-nufft_multi-coils.svg", format="svg", height=400, width=800)

