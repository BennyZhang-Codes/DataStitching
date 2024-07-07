using KomaHighOrder
seq = read_seq("E:/skope/skope_pulseq/xw_sp2d-1mm-r1_noDUM.seq") # skope sequence
grad = MAT.matread("E:/skope/20240308/grad_1mm.mat") # skope measured gradients
Δt = grad["dt"]
skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3  # mT
skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3  # mT
 
# t = range(0, 88.1*1e-3; length=88101)
t = Δt * ones(88100)
GR_skope = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1)

# hoseq = HO_Sequence(seq[8], GR_skope, maximum(GR_skope.dur))
# hoseqd = HO_discretize(hoseq)
hoseq = HO_Sequence(seq)
hoseq.GR_skope[:,8] = GR_skope
HO_plot_hoseqd(hoseq)
HO_plot_seq(hoseq)


sim_params=KomaMRICore.default_sim_params()
t, Δt      = KomaMRIBase.get_variable_times(seq; Δt=sim_params["Δt"], Δt_rf=sim_params["Δt_rf"])
B1, Δf     = KomaMRIBase.get_rfs(seq, t)
Gx, Gy, Gz = KomaMRIBase.get_grads(seq, t)
tadc       = KomaMRIBase.get_adc_sampling_times(seq)
ADCflag    = [any(tt .== tadc) for tt in t]
seqd = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)  


