using KomaHighOrder
seq = read_seq("E:/skope/skope_pulseq/xw_sp2d-1mm-r1_noDUM.seq") # skope sequence
grad = MAT.matread("E:/skope/20240308/grad_1mm.mat") # skope measured gradients
Δt = grad["dt"]
skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3  # mT
skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3  # mT
 
# t = range(0, 88.1*1e-3; length=88101)
t = Δt * ones(88100)
GR_skope = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1)

hoseq = HO_Sequence(seq)
hoseq.GR_skope[:,8] = GR_skope
# hoseqd = discretize(hoseq)

obj = brain_phantom2D(brain3D_02(),0.8;ss=10)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
raw = simulate(obj, hoseq, sys)



# simulate
obj = brain_phantom2D(ss=4)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"] = BlochHighOrder()
sim_params["Nblocks"] = 50
raw = simulate(obj, hoseq, sys; sim_params)



