################################
# with coil sensitivity
################################

# 1. hoseq
hoseq = demo_hoseq();
# plot_hoseqd(hoseq);

# 2. phantom
obj = brain_phantom2D(brain2D(); ss=10, location=0.8, coil_type=:birdcage, Nparts=5); info(obj)
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π  cancel Δw
plot_phantom_map(obj, :Dθ; darkmode=true)

# ref = brain_phantom2D_reference(brain2D(); B0map=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
# p_Δw_ref = plot_image(ref; title="quadraticB0map, [-$maxOffresonance,$maxOffresonance] Hz", darkmode=true, zmin=-maxOffresonance)


# 3. scanner & sim_params
sys = Scanner();
sim_method=BlochHighOrder("000");
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = sim_method;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal);
image = reconstruct_2d_image(raw);
p_image = plot_image(image; darkmode=true, title="Sim: $(sim_method.name)")
# savefig(p_image, (@__DIR__)*"/Sim000_quadraticB0map$maxOffresonance.svg", width=550,height=500,format="svg")




# 1. hoseq
hoseq = demo_hoseq();
# plot_hoseqd(hoseq);

# 2. phantom
maxOffresonance = 5.
obj = brain_phantom2D(brain2D(); ss=3, location=0.8, B0map=:quadratic, maxOffresonance=maxOffresonance); info(obj)
# obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π  cancel Δw
# obj.T2 .= obj.T2 * Inf;   # cancel T2 relaxiation
p_Δw = plot_phantom_map(obj, :Δw; darkmode=true)
# savefig(p_Δw, (@__DIR__)*"/quadraticB0map_Δw_$maxOffresonance.svg", width=500,height=500,format="svg")
ref = brain_phantom2D_reference(brain2D(); B0map=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
p_Δw_ref = plot_image(ref; title="quadraticB0map, [-$maxOffresonance,$maxOffresonance] Hz", darkmode=true, zmin=-maxOffresonance)
# savefig(p_Δw_ref, (@__DIR__)*"/quadraticB0map_$maxOffresonance.svg", width=550,height=500,format="svg")

# 3. scanner & sim_params
sys = Scanner();
sim_method::BlochHighOrder=BlochHighOrder("000");
sim_params = KomaMRICore.default_sim_params(sim_params)
sim_params["sim_method"] = sim_method;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal);
image = reconstruct_2d_image(raw);
p_image = plot_image(image; darkmode=true, title="Sim: 000, Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, (@__DIR__)*"/Sim000_quadraticB0map$maxOffresonance.svg", width=550,height=500,format="svg")



# 1. hoseq
hoseq = demo_hoseq();
# plot_hoseqd(hoseq);

# 2. phantom
# maxOffresonance = 20.
obj = brain_phantom2D(brain2D(); ss=3, location=0.8, B0map=:file, maxOffresonance=maxOffresonance); info(obj)
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π  cancel Δw
# obj.T2 .= obj.T2 * Inf;   # cancel T2 relaxiation
p_Δw = plot_phantom_map(obj, :Δw; darkmode=true)
# savefig(p_Δw, (@__DIR__)*"/brain2D_B0map_Δw.svg", width=500,height=500,format="svg")
ref = brain_phantom2D_reference(brain2D(); B0map=:file,key=:Δw, maxOffresonance=maxOffresonance); 
p_Δw_ref = plot_image(ref; title="brain2D_B0map, [$(round(minimum(ref))),$(round(maximum(ref)))] Hz", darkmode=true, zmin=-270)
# savefig(p_Δw_ref, (@__DIR__)*"/brain2D_B0map.svg", width=550,height=500,format="svg")

# 3. scanner & sim_params
sys = Scanner();
sim_method::BlochHighOrder=BlochHighOrder("000");
sim_params = KomaMRICore.default_sim_params(sim_params)
sim_params["sim_method"] = sim_method;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal);
image = reconstruct_2d_image(raw);
p_image = plot_image(image; darkmode=true, title="Sim: 000, Δw: [$(round(minimum(ref))),$(round(maximum(ref)))] Hz")
# savefig(p_image, (@__DIR__)*"/Sim000_brain2D_B0map.svg", width=550,height=500,format="svg")



