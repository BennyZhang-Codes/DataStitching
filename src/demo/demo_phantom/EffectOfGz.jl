using CUDA
device!(1) 
path = @__DIR__
using KomaHighOrder


obj = brain_phantom3D(brain3D(file = phantom_dict[:brain3d_285]);ss=1, start_end=[1,400]); info(obj);
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π

BHO_name = "000";
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(BHO_name);
sim_params["Nblocks"] = 3000;
sim_params["gpu"] = true;
sim_params["gpu_device"] = 0;
sim_params["return_type"]="mat";

########################################################################################################################
# without h3
########################################################################################################################
seq = demo_seq()
hoseq = HO_Sequence(seq);
# plot_seq(hoseq)

signal = simulate(obj, hoseq, sys; sim_params);
raw_withoutGz = signal_to_raw_data(signal, hoseq, :nominal);

plot_image(reconstruct_2d_image(raw_withoutGz); title="$(sim_params["sim_method"]) withoutGz Nominal", height=700, width=750)
protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal"
raw_withoutGz.params["protocolName"] = protocolName
mrd = ISMRMRDFile(path*"/src/demo/demo_phantom/rawdata/$(protocolName)_brain3d_285_withoutGz.mrd")
save(mrd, raw_withoutGz)



########################################################################################################################
# with h3
########################################################################################################################
seq = demo_seq()
GR_skope = demo_GR_skope();
seq.GR[3,8] = GR_skope[4,1]
hoseq = HO_Sequence(seq);
plot_seq(hoseq;darkmode=true)


signal = simulate(obj, hoseq, sys; sim_params);
raw_withGz = signal_to_raw_data(signal, hoseq, :nominal);

plot_image(reconstruct_2d_image(raw_withGz); title="$(sim_params["sim_method"]) withGz Nominal", height=700, width=750)
protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)_nominal"
raw_withGz.params["protocolName"] = protocolName
mrd = ISMRMRDFile(path*"/src/demo/demo_phantom/rawdata/$(protocolName)_brain3d_285_withGz.mrd")
save(mrd, raw_withGz)


img_withGz = reconstruct_2d_image(raw_withGz)
img_withoutGz = reconstruct_2d_image(raw_withoutGz)
img_error = img_withGz - img_withoutGz


width = 400
height = 350

p_withGz = plot_image(img_withGz; title="withGz Nominal", width=width, height=height)
p_withoutGz = plot_image(img_withoutGz; title="withoutGz Nominal", width=width, height=height)
p_error = plot_image(img_error; title="withGz - withoutGz", width=width, height=height)

savefig(p_withGz,       path*"/throughplane_withGz.svg", width=width, height=height,format="svg")
savefig(p_withoutGz,     path*"/throughplane_withoutGz.svg", width=width, height=height,format="svg")
savefig(p_error,      path*"/throughplane_difference.svg", width=width, height=height,format="svg")
