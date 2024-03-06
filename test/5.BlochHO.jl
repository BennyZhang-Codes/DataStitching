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

p = brain_phantom2D(ss=10)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()


hoseqd = discretize(hoseq)
sim_method = BlochHighOrder()
seq = hoseqd.seqd
#Simulation
#Motion
xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
#Effective field
Bzh0 = hoseqd.h0'
Bzh1 = xt .* hoseqd.h1'
Bzh2 = yt .* hoseqd.h2'
Bzh3 = zt .* hoseqd.h3'
Bzh4 = (xt .* yt) .* hoseqd.h4'
Bzh5 = (zt .* yt) .* hoseqd.h5'
Bzh6 = (3zt.^2-(xt.^2 .+ yt.^2 .+ zt.^2)) .* hoseqd.h6'
Bzh7 = (xt .* zt) .* hoseqd.h7'
Bzh8 = (xt.^2 .+ yt.^2) .* hoseqd.h8'

Bzho1 = Bzh1 .+ Bzh2 .+ Bzh3
Bzho2 = Bzh4 .+ Bzh5 .+ Bzh6 .+ Bzh7 .+ Bzh8

Bz0 = sim_method.ho0 ? Bzh0 : 0
Bz1 = sim_method.ho1 ? Bzho1 : xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz'
Bz2 = sim_method.ho2 ? Bzho2 : 0


p0 = plot(Bzh0[1,:]*1e3)
p1 = plot(Bzh1[1000,:]*1e3)
p2 = plot(Bzh2[1000,:]*1e3)
p3 = plot(Bzh3[1000,:]*1e3)
p4 = plot(Bzh4[1000,:]*1e3)
p5 = plot(Bzh5[1000,:]*1e3)
p6 = plot(Bzh6[1000,:]*1e3)
p7 = plot(Bzh7[1000,:]*1e3)
p8 = plot(Bzh8[1000,:]*1e3)

using PlotlyJS
data = [Bzh0[1,:]*1e3 Bzh1[1000,:]*1e3 Bzh2[1000,:]*1e3 Bzh3[1000,:]*1e3 Bzh4[1000,:]*1e3 Bzh5[1000,:]*1e3 Bzh6[1000,:]*1e3 Bzh7[1000,:]*1e3 Bzh8[1000,:]*1e3]
plot(data)


plot([(xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz')[1000,:] Bzho1[1000,:]]*1e3)
