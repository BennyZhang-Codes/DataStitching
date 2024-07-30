using KomaHighOrder
path = @__DIR__
path = path*"/src/demo/"
seq = read_seq(path*"/xw_sp2d-1mm-r1_noDUM.seq"); # skope sequence
grad = MAT.matread(path*"/grad_1mm.mat"); # skope measured gradients

# skope 
Δt = grad["dt"];
skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3; 
skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3;
t = Δt * ones(88100);
GR_dfc = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);

# hoseq
seq.GR[1,:] = -seq.GR[1,:];
hoseq = HO_Sequence(seq);
hoseq.GR_dfc[2:4, :] = hoseq.SEQ.GR;
hoseq.GR_dfc[:,8] = GR_dfc;
plot_hoseqd(hoseq)

# phantom
obj = brain_phantom2D(brain3D_02(); ss=3, location=0.8); info(obj);
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π


# scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(ho0=false, ho1=true, ho2=false);
sim_params["Nblocks"] = 150;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

# simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw_nominal = signal_to_raw_data(signal, hoseq, :nominal)
raw_measured = signal_to_raw_data(signal, hoseq, :measured)

plot_image(reconstruct_2d_image(raw_nominal); title="$(sim_params["sim_method"]) Nominal", height=700, width=750)
plot_image(reconstruct_2d_image(raw_measured); title="$(sim_params["sim_method"]) Measured", height=700, width=750)

using Dates
mrd = Sys.iswindows() ? ISMRMRDFile("E:/skope/$(seq.DEF["Name"])_$(Dates.now()).mrd") : ISMRMRDFile("/home/jyzhang/Desktop/skope/$(seq.DEF["Name"])_$(Dates.now()).mrd")
save(mrd, raw)


# gre reference
seq = read_seq("/home/jyzhang/Desktop/skope/20240104_gre/DEMO_gre2_TR1000_TE50_FA54_150_FOV150.seq");
seq.GR[1,:] = -seq.GR[1,:];
hoseq = HO_Sequence(seq);
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(ho0=false, ho1=false, ho2=false);
sim_params["Nblocks"] = 150;
sim_params["gpu"] = true;

raw = simulate(obj, hoseq, sys; sim_params);
plot_image(reconstruct_2d_image(raw); title="$(sim_params["sim_method"])", height=700, width=750)


seq = read_seq("/home/jyzhang/Desktop/skope/20240202_spiral/spiral_fov220_nx220_1.00.seq");
seq.GR[1,:] = -seq.GR[1,:];
hoseq = HO_Sequence(seq);
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(ho0=false, ho1=false, ho2=false);
sim_params["Nblocks"] = 150;
sim_params["gpu"] = true;

raw = simulate(obj, hoseq, sys; sim_params);
plot_image(reconstruct_2d_image(raw); title="$(sim_params["sim_method"])", height=700, width=750)

