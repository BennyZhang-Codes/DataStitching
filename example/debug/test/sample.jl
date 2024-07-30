using KomaHighOrder
if Sys.iswindows()
    seq = read_seq("E:/skope/skope_pulseq/xw_sp2d-1mm-r1_noDUM.seq"); # skope sequence
    grad = MAT.matread("E:/skope/20240308/grad_1mm.mat"); # skope measured gradients
elseif Sys.islinux()
    seq = read_seq("/home/jyzhang/Desktop/skope/20240308/xw_sp2d-1mm-r1_noDUM.seq"); # skope sequence
    grad = MAT.matread("/home/jyzhang/Desktop/skope/20240308/grad_1mm.mat"); # skope measured gradients
end
# skope 
Δt = grad["dt"];
skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3; 
skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3;

t = Δt * ones(88100);
GR_dfc = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);

# hoseq
seq.GR[1,:] = -seq.GR[1,:];
hoseq = HO_Sequence(seq);
hoseq.GR_dfc[:,8] = GR_dfc;

obj = brain_phantom2D(brain3D_02(); ss=10, location=0.8); info(obj);
obj.Δw .= 0;

sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(ho0=false, ho1=false, ho2=false);
sim_params["Nblocks"] = 150;
sim_params["gpu"] = true;
raw = simulate(obj, hoseq, sys; sim_params);

# using Dates
# mrd = Sys.iswindows() ? ISMRMRDFile("E:/skope/$(seq.DEF["Name"])_$(Dates.now()).mrd") : ISMRMRDFile("/home/jyzhang/Desktop/skope/$(seq.DEF["Name"])_$(Dates.now()).mrd")
# save(mrd, raw)

# plot_signal(raw)
img = recon_2d(raw)
p = plot_image(img; title="$(sim_params["sim_method"])", height=600, width=750)

plot_phantom_map(obj, :Δw; darkmode=true, width=500, height=500)

