import KomaMRI.KomaMRICore: simulate
"""
    out = simulate(obj::HO_Phantom, hoseq::HO_Sequence, sys::Scanner; sim_params, w)

Returns the raw signal or the last state of the magnetization according to the value
of the `"return_type"` key of the `sim_params` dictionary.

# Arguments
- `obj`: (`::HO_Phantom`) HO_Phantom struct
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `sim_params`: (`::Dict{String,Any}`, `=Dict{String,Any}()`) simulation parameter dictionary
- `w`: (`::Blink.AtomShell.Window`, `=nothing`) the window within which to display a
    progress bar in the Blink Window UI. If this variable is anything other than 'nothing',
    the progress bar will be considered

# Returns
- `out`: (`::Vector{Complex}` or `::SpinStateRepresentation` or `::RawAcquisitionData`) depending
    on whether "return_type" is "mat", "state" or "raw" (default), respectively
"""
function simulate(
    obj::HO_Phantom, hoseq::HO_Sequence, sys::Scanner;
    sim_params=Dict{String,Any}(), w=nothing
) 
    #Simulation parameter unpacking, and setting defaults if key is not defined
    sim_params = KomaMRICore.default_sim_params(sim_params)
    # Simulation init
    hoseqd = discretize(hoseq; sampling_params=sim_params) # Sampling of Sequence waveforms
    parts, excitation_bool = KomaMRICore.get_sim_ranges(hoseqd.seqd; Nblocks=sim_params["Nblocks"]) # Generating simulation blocks
    t_sim_parts = [hoseqd.seqd.t[p[1]] for p ∈ parts]; append!(t_sim_parts, hoseqd.seqd.t[end])
    # Spins' state init (Magnetization, EPG, etc.), could include modifications to obj (e.g. T2*)
    Xt, obj = initialize_spins_state(obj, sim_params["sim_method"])
    # Signal init
    Ndims = sim_output_dim(obj, hoseq.SEQ, sys, sim_params["sim_method"])
    sig = zeros(ComplexF64, Ndims..., sim_params["Nthreads"])
    # Objects to GPU
    if sim_params["gpu"] #Default
        device!(sim_params["gpu_device"])
        gpu_name = name.(devices())[sim_params["gpu_device"]+1]
        obj    = obj    |> gpu #HO_Phantom
        hoseqd = hoseqd |> gpu #HO_DiscreteSequence
        Xt     = Xt     |> gpu #SpinStateRepresentation
        sig    = sig    |> gpu #Signal
    end
    if sim_params["precision"] == "f32" #Default
        obj    = obj    |> f32 #HO_Phantom
        hoseqd = hoseqd |> f32 #HO_DiscreteSequence
        Xt     = Xt     |> f32 #SpinStateRepresentation
        sig    = sig    |> f32 #Signal
    elseif sim_params["precision"] == "f64"
        obj    = obj    |> f64 #HO_Phantom
        hoseqd = hoseqd |> f64 #HO_DiscreteSequence
        Xt     = Xt     |> f64 #SpinStateRepresentation
        sig    = sig    |> f64 #Signal
    end
    # Simulation
    @info "Running simulation in the $(sim_params["gpu"] ? "GPU ($gpu_name)" : "CPU with $(sim_params["Nthreads"]) thread(s)")" koma_version=KomaMRICore.__VERSION__ sim_method = sim_params["sim_method"] spins = length(obj) time_points = length(hoseqd.seqd.t) adc_points=Ndims[1]
    @time timed_tuple = @timed run_sim_time_iter!(obj, hoseqd, sig, Xt, sim_params["sim_method"]; Nblocks=length(parts), Nthreads=sim_params["Nthreads"], parts, excitation_bool, w)
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
        out = signal_to_raw_data(sig, hoseq; phantom_name=obj.name, sys=sys, sim_params=sim_params_raw)
    end
    return out
end