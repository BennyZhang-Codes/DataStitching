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
# hoseqd = HO_discretize(hoseq)

obj = brain_phantom2D(ss=10)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
raw = HO_simulate(obj, hoseq, sys)


using CUDA
# Printing
using ProgressMeter

obj = brain_phantom2D(ss=10)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
w=nothing

#Simulation parameter unpacking, and setting defaults if key is not defined
sim_params = KomaMRICore.default_sim_params(sim_params)
# Simulation init
hoseqd = HO_discretize(hoseq; sampling_params=sim_params) # Sampling of Sequence waveforms
parts, excitation_bool = KomaMRICore.get_sim_ranges(hoseqd.seqd; Nblocks=sim_params["Nblocks"]) # Generating simulation blocks
t_sim_parts = [hoseqd.seqd.t[p[1]] for p ∈ parts]; append!(t_sim_parts, hoseqd.seqd.t[end])
# Spins' state init (Magnetization, EPG, etc.), could include modifications to obj (e.g. T2*)
Xt, obj = KomaMRICore.initialize_spins_state(obj, sim_params["sim_method"])
# Signal init
Ndims = KomaMRICore.sim_output_dim(obj, hoseq.SEQ, sys, sim_params["sim_method"])
sig = zeros(ComplexF64, Ndims..., sim_params["Nthreads"])
# Objects to GPU
if sim_params["gpu"] #Default
    device!(sim_params["gpu_device"])
    gpu_name = name.(devices())[sim_params["gpu_device"]+1]
    obj    = obj    |> gpu #Phantom
    hoseqd = hoseqd |> gpu #DiscreteSequence
    Xt     = Xt     |> gpu #SpinStateRepresentation
    sig    = sig    |> gpu #Signal
end
if sim_params["precision"] == "f32" #Default
    obj    = obj    |> f32 #Phantom
    hoseqd = hoseqd |> f32 #DiscreteSequence
    Xt     = Xt     |> f32 #SpinStateRepresentation
    sig    = sig    |> f32 #Signal
elseif sim_params["precision"] == "f64"
    obj    = obj    |> f64 #Phantom
    hoseqd = hoseqd |> f64 #DiscreteSequence
    Xt     = Xt     |> f64 #SpinStateRepresentation
    sig    = sig    |> f64 #Signal
end
# Simulation
@info "Running simulation in the $(sim_params["gpu"] ? "GPU ($gpu_name)" : "CPU with $(sim_params["Nthreads"]) thread(s)")" koma_version=KomaMRICore.__VERSION__ sim_method = sim_params["sim_method"] spins = length(obj) time_points = length(hoseqd.seqd.t) adc_points=Ndims[1]
@time timed_tuple = @timed KomaHighOrder.run_sim_time_iter!(obj, hoseqd, sig, Xt, sim_params["sim_method"]; Nblocks=length(parts), Nthreads=sim_params["Nthreads"], parts, excitation_bool, w)
# Result to CPU, if already in the CPU it does nothing
sig = sum(sig; dims=length(Ndims)+1) |> cpu
sig .*= KomaMRICore.get_adc_phase_compensation(hoseq.SEQ)
Xt = Xt |> cpu
if sim_params["gpu"] GC.gc(true); CUDA.reclaim() end
# Output
if sim_params["return_type"] == "state"
    out = Xt
elseif sim_params["return_type"] == "mat"
    out = sig
elseif sim_params["return_type"] == "raw"
    # To visually check the simulation blocks
    sim_params_raw = copy(sim_params)
    sim_params_raw["sim_method"] = string(sim_params["sim_method"])
    sim_params_raw["gpu"] = sim_params["gpu"]
    sim_params_raw["Nthreads"] = sim_params["Nthreads"]
    sim_params_raw["t_sim_parts"] = t_sim_parts
    sim_params_raw["type_sim_parts"] = excitation_bool
    sim_params_raw["Nblocks"] = length(parts)
    sim_params_raw["sim_time_sec"] = timed_tuple.time
    out = HO_signal_to_raw_data(sig, hoseq.SEQ; phantom_name=obj.name, sys=sys, sim_params=sim_params_raw)
end


# simulate
obj = brain_phantom2D(ss=1)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
sim_params["Nblocks"] = 50
raw = HO_simulate(obj, hoseq, sys; sim_params)



