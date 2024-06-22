using KomaHighOrder

simtype  = SimType(B0=false, T2=false, ss=5)
coil_type= :birdcage
overlap  = 0
Nparts   = 9; nrows=3; ncols=3;
Npartsx  = 3
Npartsy  = 4
Ncoils = coil_type == :rect ? Npartsx * Npartsy : Nparts
BHO_name = "000"
if coil_type == :fan
    folder   = "$(coil_type)_Ncoils$(Ncoils)_overlap$(overlap)"
elseif coil_type == :rect
    folder   = "$(coil_type)_$(Npartsx)_$(Npartsy)"
elseif coil_type == :birdcage
    folder   = "$(coil_type)_Ncoils$(Ncoils)"
end
path     = "$(@__DIR__)/src/demo/demo_sense/$folder"
if ispath(path) == false mkdir(path) end

hoseq = demo_hoseq()

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

coil_images = Array{Float32, 3}(undef, Nx, Ny, Ncoils);
signal = zeros(ComplexF64, sum(hoseq.SEQ.ADC.N), Ncoils);
for coil_idx = 1:Ncoils
    obj = brain_phantom2D(brain2D(); ss=simtype.ss, location=0.8, coil_type=coil_type, coil_idx=coil_idx, Nparts=Ncoils, Npartsx=Npartsx, Npartsy=Npartsy, overlap=overlap); 
    obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; 

    # simulate
    signal[:, coil_idx] = simulate(obj, hoseq, sys; sim_params);
    protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_$(Ncoils)-$(coil_idx)"
    coil_images[:,:,coil_idx] = reconstruct_2d_image(signal_to_raw_data(signal[:, coil_idx:coil_idx], hoseq, :nominal));

    p = plot_image(coil_images[:,:,coil_idx]; title="$(protocolName)", height=400, width=450)
    savefig(p,  "$(path)/$(protocolName).svg",format="svg", height=400, width=450)
end

p_coil_images = plot_imgs_subplots(abs.(coil_images), nrows, ncols; title="$(Ncoils) coils: NUFFT recon", height=400, width=450)
savefig(p_coil_images,  "$(path)/$(raw.params["protocolName"])-nufft_multi-coils.svg", format="svg", height=400, width=450)

raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_Ncoils$(Ncoils)"
raw.params["protocolName"] = filename
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)


