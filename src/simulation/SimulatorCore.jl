
"""
    sig, Xt = run_spin_precession_parallel(obj, hoseqd, M; Nthreads)

Implementation in multiple threads for the simulation in free precession,
separating the spins of the phantom `obj` in `Nthreads`.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct

# Keywords
- `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
    dividing the simulation into different phantom spin parts

# Returns
- `sig`: (`Vector{ComplexF64}`) raw signal over time
- `Xt`: (`::Vector{Mag}`) final state of the Mag vector (or the initial state for the
    next simulation step (the next step can be another precession step or an excitation
    step))
"""
function run_spin_precession_parallel!(obj::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nthreads=Threads.nthreads()) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads)
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

    ThreadsX.foreach(enumerate(parts)) do (i, p)
        run_spin_precession!(@view(obj[p]), hoseqd, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
    end

    return nothing
end

"""
    M0 = run_spin_excitation_parallel(obj, hoseqd, Xt; Nthreads)

It gives rise to a rotation of M0 with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis. It uses
different number threads to excecute the process.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct

# Keywords
- `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
    dividing the simulation into different phantom spin parts

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (or the initial
    state for the next precession simulation step)
"""
function run_spin_excitation_parallel!(obj::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nthreads=Threads.nthreads()) where {T<:Real}

    parts = kfoldperm(length(obj), Nthreads)
    dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times

    ThreadsX.foreach(enumerate(parts)) do (i, p)
        run_spin_excitation!(@view(obj[p]), hoseqd, @view(sig[dims...,i]), @view(Xt[p]), sim_method)
    end

    return nothing
end


"""
    S_interp, M0 = run_sim_time_iter(obj, hoseqd, t, Δt; Nblocks, Nthreads, gpu, w)

Performs the simulation over the total time vector `t` by dividing the time into `Nblocks`
parts to reduce RAM usage and spliting the spins of the phantom `obj` into `Nthreads` to
take advantage of CPU parallel processing.

# Arguments
- `obj`: (`::Phantom`) Phantom struct
- `seq`: (`::Sequence`) Sequence struct
- `t`: (`::Vector{Float64}`, `[s]`) non-uniform time vector
- `Δt`: (`::Vector{Float64}`, `[s]`) delta time of `t`

# Keywords
- `Nblocks`: (`::Int`, `=16`) number of groups for spliting the simulation over time
- `Nthreads`: (`::Int`, `=Threads.nthreads()`) number of process threads for
    dividing the simulation into different phantom spin parts
- `gpu`: (`::Function`) function that represents the gpu of the host
- `w`: (`::Any`, `=nothing`) flag to regard a progress bar in the blink window UI. If
    this variable is differnet from nothing, then the progress bar is considered

# Returns
- `S_interp`: (`::Vector{ComplexF64}`) interpolated raw signal
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_sim_time_iter!(obj::Phantom, hoseqd::HO_DiscreteSequence, sig::AbstractArray{Complex{T}},
    Xt::SpinStateRepresentation{T}, sim_method::SimulationMethod;
    Nblocks=1, Nthreads=Threads.nthreads(), parts=[1:length(hoseqd)], excitation_bool=ones(Bool, size(parts)), w=nothing) where {T<:Real}
    # Simulation
    rfs = 0
    samples = 1
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        seq_block = @view hoseqd[p]
        # Params
        # excitation_bool = is_RF_on(seq_block) #&& is_ADC_off(seq_block) #PATCH: the ADC part should not be necessary, but sometimes 1 sample is identified as RF in an ADC block
        Nadc = sum(seq_block.seqd.ADC)
        acq_samples = samples:samples+Nadc-1
        dims = [Colon() for i=1:output_Ndim(sim_method)] # :,:,:,... Ndim times
        # Simulation wrappers
        if excitation_bool[block]
            run_spin_excitation_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method; Nthreads)
            rfs += 1
        else
            run_spin_precession_parallel!(obj, seq_block, @view(sig[acq_samples, dims...]), Xt, sim_method; Nthreads)
        end
        samples += Nadc
        #Update progress
        next!(progress_bar, showvalues=[(:simulated_blocks, block), (:rf_blocks, rfs), (:acq_samples, samples-1)])
        KomaMRICore.update_blink_window_progress!(w, block, Nblocks)
    end
    return nothing
end