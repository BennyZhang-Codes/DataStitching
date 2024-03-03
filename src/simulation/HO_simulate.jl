"""
    out = HO_simulate(obj::Phantom, seq::Sequence, sys::Scanner; sim_params, w)

Returns the raw signal or the last state of the magnetization according to the value
of the `"return_type"` key of the `sim_params` dictionary.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct
- `sys`: (`::Scanner`) Scanner struct

# Keywords
- `sim_params`: (`::Dict{String,Any}`, `=Dict{String,Any}()`) simulation parameter dictionary
- `w`: (`::Blink.AtomShell.Window`, `=nothing`) the window within which to display a
    progress bar in the Blink Window UI. If this variable is anything other than 'nothing',
    the progress bar will be considered

# Returns
- `out`: (`::Vector{Complex}` or `::SpinStateRepresentation` or `::RawAcquisitionData`) depending
    on whether "return_type" is "mat", "state" or "raw" (default), respectively

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/3.koma_paper/comparison_accuracy/sequences/EPI/epi_100x100_TE100_FOV230.seq");

julia> sys, obj, seq = Scanner(), brain_phantom2D(), read_seq(seq_file)

julia> raw = simulate(obj, seq, sys)

julia> plot_signal(raw)
```
"""
function HO_simulate(
    obj::Phantom, seq::Sequence, sys::Scanner;
    sim_params=Dict{String,Any}(), w=nothing
)
    #Simulation parameter unpacking, and setting defaults if key is not defined
    sim_params = KomaMRICore.default_sim_params(sim_params)
    # Simulation init
    seqd = KomaMRICore.discretize(seq; sampling_params=sim_params) # Sampling of Sequence waveforms
    parts, excitation_bool = KomaMRICore.get_sim_ranges(seqd; Nblocks=sim_params["Nblocks"]) # Generating simulation blocks
    t_sim_parts = [seqd.t[p[1]] for p ∈ parts]; append!(t_sim_parts, seqd.t[end])
    # Spins' state init (Magnetization, EPG, etc.), could include modifications to obj (e.g. T2*)
    Xt, obj = KomaMRICore.initialize_spins_state(obj, sim_params["sim_method"])
    # Signal init
    Ndims = KomaMRICore.sim_output_dim(obj, seq, sys, sim_params["sim_method"])
    sig = zeros(ComplexF64, Ndims..., sim_params["Nthreads"])
    # Objects to GPU
    if sim_params["gpu"] #Default
        device!(sim_params["gpu_device"])
        gpu_name = name.(devices())[sim_params["gpu_device"]+1]
        obj  = obj  |> gpu #Phantom
        seqd = seqd |> gpu #DiscreteSequence
        Xt   = Xt   |> gpu #SpinStateRepresentation
        sig  = sig  |> gpu #Signal
    end
    if sim_params["precision"] == "f32" #Default
        obj  = obj  |> f32 #Phantom
        seqd = seqd |> f32 #DiscreteSequence
        Xt   = Xt   |> f32 #SpinStateRepresentation
        sig  = sig  |> f32 #Signal
    elseif sim_params["precision"] == "f64"
        obj  = obj  |> f64 #Phantom
        seqd = seqd |> f64 #DiscreteSequence
        Xt   = Xt   |> f64 #SpinStateRepresentation
        sig  = sig  |> f64 #Signal
    end
    # Simulation
    @info "Running simulation in the $(sim_params["gpu"] ? "GPU ($gpu_name)" : "CPU with $(sim_params["Nthreads"]) thread(s)")" koma_version=__VERSION__ sim_method = sim_params["sim_method"] spins = length(obj) time_points = length(seqd.t) adc_points=Ndims[1]
    @time timed_tuple = @timed KomaMRICore.run_sim_time_iter!(obj, seqd, sig, Xt, sim_params["sim_method"]; Nblocks=length(parts), Nthreads=sim_params["Nthreads"], parts, excitation_bool, w)
    # Result to CPU, if already in the CPU it does nothing
    sig = sum(sig; dims=length(Ndims)+1) |> cpu
    sig .*= KomaMRICore.get_adc_phase_compensation(seq)
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
        out = HO_signal_to_raw_data(sig, seq; phantom_name=obj.name, sys=sys, sim_params=sim_params_raw)
    end
    return out
end