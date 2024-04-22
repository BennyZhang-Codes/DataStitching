using KomaHighOrder


Ncoils = 30
BHO_name = "000"


folder = "Ncoils$Ncoils"
path = "$(@__DIR__)/src/demo/demo_sense/$folder"
if ispath(path) == false mkdir(path) end

hoseq = demo_hoseq()

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["gpu"] = true;
sim_params["return_type"]="mat";


signal = zeros(ComplexF64, sum(hoseq.SEQ.ADC.N), Ncoils);
for coil_idx = 1:Ncoils
    obj = brain_phantom2D(brain2D(); ss=3, location=0.8, mask_idx=coil_idx, Nparts=Ncoils); 
    obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    # obj.T2 .= obj.T2 * Inf; 
    # simulate
    signal[:, coil_idx] = simulate(obj, hoseq, sys; sim_params);
    protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_$(Ncoils)-$(coil_idx)"
    # p = plot_image(reconstruct_2d_image(signal_to_raw_data(signal[:, coil_idx:coil_idx], hoseq, :nominal)); 
    #     title="$(protocolName)", height=400, width=450)
    # savefig(p,  "$(path)/$(protocolName).svg",format="svg", height=400, width=450)
end

raw = signal_to_raw_data(signal, hoseq, :nominal)
filename = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal_Ncoils$(Ncoils)"
raw.params["protocolName"] = filename
mrd = ISMRMRDFile("$(path)/$(filename).mrd")
save(mrd, raw)


