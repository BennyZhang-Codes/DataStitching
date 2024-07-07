using KomaHighOrder, PlotlyJS

hoseq = demo_hoseq()
p = brain_phantom2D(BrainPhantom(), ss=10, csmtype=:fan, nCoil=3, overlap=0)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()


hoseqd = discretize(hoseq)[80000:85000]
seq = hoseqd.seqd
sim_method = BlochHighOrder()
M, _ = KomaMRICore.initialize_spins_state(p, sim_method)
sig = zeros(ComplexF64, sum(seq.ADC))
#Simulation
#Coil sensitivity
smap = p.Dλ1 .+ p.Dλ2*im
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

Bz = Bz0 .+ Bz1 .+ Bz2 .+ p.Δw / (2π * γ)
#Rotation
if is_ADC_on(seq)
    ϕ = (-2π * γ) .* cumtrapz(seq.Δt', Bz)
else
    ϕ = (-2π * γ) .* trapz(seq.Δt', Bz)
end
#Mxy precession and relaxation, and Mz relaxation
tp = cumsum(seq.Δt) # t' = t - t0
dur = sum(seq.Δt)   # Total length, used for signal relaxation
Mxy = [M.xy M.xy .* exp.(1im .* ϕ .- tp' ./ p.T2)] #This assumes Δw and T2 are constant in time
M.xy .= Mxy[:, end]
M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
sig .= transpose(sum(Mxy[:, findall(seq.ADC)] .* smap; dims=1))


using PlotlyJS
data = [Bzh0[1,:]*1e3 Bzh1[1000,:]*1e3 Bzh2[1000,:]*1e3 Bzh3[1000,:]*1e3 Bzh4[1000,:]*1e3 Bzh5[1000,:]*1e3 Bzh6[1000,:]*1e3 Bzh7[1000,:]*1e3 Bzh8[1000,:]*1e3]
plot(data)


plot([(xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz')[1000,:] Bzho1[1000,:]]*1e3)
