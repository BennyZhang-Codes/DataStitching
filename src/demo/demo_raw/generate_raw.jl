using KomaHighOrder
hoseq = demo()
plot_hoseqd(hoseq)

# phantom
obj = brain_phantom2D(brain3D_02(); ss=3, location=0.8); info(obj);
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π

BHO_name = "000"
# scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
# sim_params["Nblocks"] = 150;
sim_params["gpu"] = true;
sim_params["gpu_device"] = 1;
sim_params["return_type"]="mat";

# simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal)
plot_image(reconstruct_2d_image(raw); title="$(sim_params["sim_method"]) Nominal", height=700, width=750)

protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal"
raw.params["protocolName"] = protocolName
path = @__DIR__
mrd = ISMRMRDFile(path*"/src/demo/demo_raw/$(protocolName).mrd")
save(mrd, raw)
