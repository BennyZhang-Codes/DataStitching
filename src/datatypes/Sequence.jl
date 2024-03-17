import KomaMRI.KomaMRIBase: get_samples, get_theo_Gi, get_grads, get_kspace, get_theo_t, get_theo_A

"""
    hoseq = HO_Sequence()
    hoseq = HO_Sequence(GR_skope::Array{Grad,2})
    hoseq = HO_Sequence(SEQ::Sequence, GR_skope::Array{Grad,2})
	hoseq = HO_Sequence(SEQ::Sequence)
	hoseq = HO_Sequence(GR_skope::Array{Grad,2}, SEQ::Sequence)
	hoseq = HO_Sequence(GR_skope::Array{Grad,1})
	hoseq = HO_Sequence(SEQ::Sequence, GR_skope::Array{Grad,1}, DUR)


The HO_Sequence struct. It contains MRI sequence and corresponding skope measured gradients.

# Arguments
- `SEQ`: (`::Sequence`) MRI sequence.
- `GR_skope`: (`::Matrix{Grad}`) gradient matrix for skope. Rows for amplitudes and columns are for blocks
\n 	h0 := 1 
\n  h1 := x
\n  h2 := y 
\n  h3 := z
\n  h4 := xy
\n  h5 := yz
\n  h6 := 3z² - (x²	+ y² + z²)
\n  h7 := xz
\n  h8 := x² - y²
- `DUR`: (`::Vector`, `[s]`) duration block vector

# Returns
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct
"""
mutable struct HO_Sequence
	SEQ::Sequence
	GR_skope::Array{Grad,2}	  #Sequence in (h0, ..., h8) and time
	DUR::Vector				  #Duration of each block, this enables delays after RF pulses to satisfy ring-down times
	HO_Sequence(SEQ, GR_skope, DUR) = begin
		@assert size(SEQ)[1] .== size(GR_skope, 2)  "The number of SEQ and GR_skope objects must be the same."
		# @assert size(GR_skope, 1) .== 9 "The number of rows in the GR_skope matrix must be 9."
		M, N = size(GR_skope)
		new(SEQ,
		    [i <= M ? GR_skope[i,j] : Grad(0, 0) for i in 1:9, j in 1:N],
			maximum([SEQ.DUR GR_skope.dur DUR],dims=2)[:])
	end
end

# Main Constructors
function HO_Sequence(GR_skope::Array{Grad,2})
	M, N = size(GR_skope)
    gr_skope = [i <= M ? GR_skope[i,j] : Grad(0, 0) for i in 1:9, j in 1:N]
	M, N = size(gr_skope[2:3,:])   
	seq = Sequence([i <= M ? gr_skope[i,j] : Grad(0, 0) for i in 1:3, j in 1:N])  # use skope gradient to generate seq
	dur = maximum([seq.DUR gr_skope.dur],dims=2)[:]
    return HO_Sequence(seq, gr_skope, dur)
end

function HO_Sequence(SEQ::Sequence, GR_skope::Array{Grad,2})
	@assert size(SEQ)[1] .== size(GR_skope, 2) "The number of SEQ and GR_skope objects must be the same."
	M, N = size(GR_skope)
    gr_skope = [i <= M ? GR_skope[i,j] : Grad(0, 0) for i in 1:9, j in 1:N]
	dur = maximum([SEQ.DUR gr_skope.dur],dims=2)[:]
    return HO_Sequence(SEQ, gr_skope, dur)
end

function HO_Sequence(SEQ::Sequence)
	gr_skope = [Grad(0, 0) for _ in 1:9, _ in 1:size(SEQ)[1]]
	dur = maximum([SEQ.DUR gr_skope.dur],dims=2)[:]
    return HO_Sequence(SEQ, gr_skope, dur)
end

# Other constructors
HO_Sequence(GR_skope::Array{Grad,1}) = HO_Sequence(reshape(GR_skope,1,:))
HO_Sequence(SEQ::Sequence, GR_skope::Array{Grad,1}, DUR)= HO_Sequence(SEQ, reshape(GR_skope, :, 1), [DUR])
HO_Sequence() = HO_Sequence([Grad(0, 0)])


#Sequence operations
Base.length(x::HO_Sequence) = length(x.DUR)
Base.iterate(x::HO_Sequence) = (HO_Sequence(x.SEQ[1], x.GR_skope[:,1], x.DUR[1]), 2)
Base.iterate(x::HO_Sequence, i::Integer) = (i <= length(x)) ? (HO_Sequence(x.SEQ[i], x.GR_skope[:,i], x.DUR[i]), i+1) : nothing
Base.getindex(x::HO_Sequence, i::UnitRange{Int}) = HO_Sequence(x.SEQ[i], x.GR_skope[:,i], x.DUR[i])
Base.getindex(x::HO_Sequence, i::Int) = HO_Sequence(x.SEQ[i], x.GR_skope[:,i], x.DUR[i])
Base.getindex(x::HO_Sequence, i::BitArray{1}) = any(i) ? HO_Sequence(x.SEQ[i], x.GR_skope[:,i], x.DUR[i]) : nothing
Base.getindex(x::HO_Sequence, i::Array{Bool,1}) = any(i) ? HO_Sequence(x.SEQ[i], x.GR_skope[:,i], x.DUR[i]) : nothing
Base.lastindex(x::HO_Sequence) = length(x.DUR)
# copy(x::HO_Sequence) where HO_Sequence = HO_Sequence(copy(x.SEQ), deepcopy(x.GR_skope), deepcopy(x.DUR))


#Arithmetic operations
+(x::HO_Sequence, y::HO_Sequence) = HO_Sequence(x.SEQ + y.SEQ, [x.GR_skope  y.GR_skope], [x.DUR; y.DUR])
-(x::HO_Sequence, y::HO_Sequence) = HO_Sequence(x.SEQ + y.SEQ, [x.GR_skope -y.GR_skope], [x.DUR; y.DUR])
-(x::HO_Sequence) = HO_Sequence(-x.SEQ, -x.GR_skope, x.DUR)
*(x::HO_Sequence, α::Real) = HO_Sequence(α*x.SEQ, α*x.GR_skope, x.DUR)
*(α::Real, x::HO_Sequence) = HO_Sequence(α*x.SEQ, α*x.GR_skope, x.DUR)
*(x::HO_Sequence, α::ComplexF64) = HO_Sequence(α*x.SEQ, x.GR_skope, x.DUR)
*(α::ComplexF64, x::HO_Sequence) = HO_Sequence(α*x.SEQ, x.GR_skope, x.DUR)
/(x::HO_Sequence, α::Real) = HO_Sequence(x.SEQ/α, x.GR_skope, x.DUR)

#Sequence object functions
size(x::HO_Sequence) = size(x.GR_skope[1,:])

"""
    str = show(io::IO, s::HO_Sequence)

Displays information about the Sequence struct `s` in the julia REPL.

# Arguments
- `s`: (`::HO_Sequence`) Sequence struct

# Returns
- `str` (`::String`) output string message
"""
Base.show(io::IO, s::HO_Sequence) = begin
	show(io, s.SEQ)
	name = ["h$(i)" for i in 0:8]
	for i in 1:9
		print(io, "\n$(name[i]): $(s.GR_skope[i, :])")
	end
end

"""
    samples = get_samples(hoseq::HO_Sequence; off_val=0, max_rf_samples=Inf)

Returns the samples of the events in `hoseq`.

# Arguments
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct

# Keywords
- `off_val`: (`::Number`, `=0`) offset value for amplitude. Typically used to hide points in
    plots by setting it to `Inf`
- `max_rf_samples`: (`::Integer`, `=Inf`) maximum number of samples for the RF struct

# Returns
- `samples`: (`::NamedTuple`) contains samples for `h0`, `h1`, `h2`, `h3`, `h4`, `h5`, `h6`,
    `h7`, `h8`, `gx`, `gy`, `gz`, `rf`, and `adc` events. Each event, represented by 
	`e::NamedTuple`, includes time samples (`e.t`) and amplitude samples (`e.A`)
"""
get_samples(hoseq::HO_Sequence; off_val=0, max_rf_samples=Inf) = begin
	gx, gy, gz, rf, adc = get_samples(hoseq.SEQ; off_val=Inf, max_rf_samples)
    N = length(hoseq.SEQ)
    T0 = get_block_start_times(hoseq.SEQ)
    # GRADs
    t_h0 = reduce(vcat, [get_theo_t(hoseq.GR_skope[1,i]) .+ T0[i] for i in 1:N])
	t_h1 = reduce(vcat, [get_theo_t(hoseq.GR_skope[2,i]) .+ T0[i] for i in 1:N])
	t_h2 = reduce(vcat, [get_theo_t(hoseq.GR_skope[3,i]) .+ T0[i] for i in 1:N])
	t_h3 = reduce(vcat, [get_theo_t(hoseq.GR_skope[4,i]) .+ T0[i] for i in 1:N])
	t_h4 = reduce(vcat, [get_theo_t(hoseq.GR_skope[5,i]) .+ T0[i] for i in 1:N])
	t_h5 = reduce(vcat, [get_theo_t(hoseq.GR_skope[6,i]) .+ T0[i] for i in 1:N])
	t_h6 = reduce(vcat, [get_theo_t(hoseq.GR_skope[7,i]) .+ T0[i] for i in 1:N])
	t_h7 = reduce(vcat, [get_theo_t(hoseq.GR_skope[8,i]) .+ T0[i] for i in 1:N])
	t_h8 = reduce(vcat, [get_theo_t(hoseq.GR_skope[9,i]) .+ T0[i] for i in 1:N])
    A_h0 = reduce(vcat, [get_theo_A(hoseq.GR_skope[1,i]; off_val) for i in 1:N])
	A_h1 = reduce(vcat, [get_theo_A(hoseq.GR_skope[2,i]; off_val) for i in 1:N])
	A_h2 = reduce(vcat, [get_theo_A(hoseq.GR_skope[3,i]; off_val) for i in 1:N])
	A_h3 = reduce(vcat, [get_theo_A(hoseq.GR_skope[4,i]; off_val) for i in 1:N])
	A_h4 = reduce(vcat, [get_theo_A(hoseq.GR_skope[5,i]; off_val) for i in 1:N])
	A_h5 = reduce(vcat, [get_theo_A(hoseq.GR_skope[6,i]; off_val) for i in 1:N])
	A_h6 = reduce(vcat, [get_theo_A(hoseq.GR_skope[7,i]; off_val) for i in 1:N])
	A_h7 = reduce(vcat, [get_theo_A(hoseq.GR_skope[8,i]; off_val) for i in 1:N])
	A_h8 = reduce(vcat, [get_theo_A(hoseq.GR_skope[9,i]; off_val) for i in 1:N])

    return (
        gx, gy, gz, rf, adc,
		h0 = (t = t_h0, A = A_h0),
		h1 = (t = t_h1, A = A_h1),
		h2 = (t = t_h2, A = A_h2),
		h3 = (t = t_h3, A = A_h3),
		h4 = (t = t_h4, A = A_h4),
		h5 = (t = t_h5, A = A_h5),
		h6 = (t = t_h6, A = A_h6),
		h7 = (t = t_h7, A = A_h7),
		h8 = (t = t_h8, A = A_h8),
    )
end

function get_theo_Gi(hoseq::HO_Sequence, idx)
    N = length(hoseq.SEQ)
    T0 = get_block_start_times(hoseq.SEQ)
    t = vcat([get_theo_t(hoseq.GR_skope[idx,i]) .+ T0[i] for i=1:N]...)
    G = vcat([get_theo_A(hoseq.GR_skope[idx,i]) for i=1:N]...) #; off_val=0 <---potential solution
	Interpolations.deduplicate_knots!(t; move_knots=true)
	return (t, G)
end


function get_grads(hoseq::HO_Sequence, t::Vector)
    h0 = get_theo_Gi(hoseq, 1)
    h1 = get_theo_Gi(hoseq, 2)
    h2 = get_theo_Gi(hoseq, 3)
    h3 = get_theo_Gi(hoseq, 4)
    h4 = get_theo_Gi(hoseq, 5)
    h5 = get_theo_Gi(hoseq, 6)
    h6 = get_theo_Gi(hoseq, 7)
    h7 = get_theo_Gi(hoseq, 8)
    h8 = get_theo_Gi(hoseq, 9)
    H0 = linear_interpolation(h0..., extrapolation_bc=0)(t)
    H1 = linear_interpolation(h1..., extrapolation_bc=0)(t)
    H2 = linear_interpolation(h2..., extrapolation_bc=0)(t)
    H3 = linear_interpolation(h3..., extrapolation_bc=0)(t)
    H4 = linear_interpolation(h4..., extrapolation_bc=0)(t)
    H5 = linear_interpolation(h5..., extrapolation_bc=0)(t)
    H6 = linear_interpolation(h6..., extrapolation_bc=0)(t)
    H7 = linear_interpolation(h7..., extrapolation_bc=0)(t)
    H8 = linear_interpolation(h8..., extrapolation_bc=0)(t)
    return (H0, H1, H2, H3, H4, H5, H6, H7, H8)
end


"""
    kspace, kspace_adc = get_kspace(seq::Sequence; Δt=1)

Outputs the designed k-space trajectory of the Sequence `seq`.

# Arguments
- `seq`: (`::Sequence`) Sequence struct
- `Δt`: (`::Real`, `=1`, `[s]`) nominal delta time separation between two time samples
    for ADC acquisition and Gradients

# Returns
- `kspace`: (`3-column ::Matrix{Float64}`) kspace
- `kspace_adc`: (`3-column ::Matrix{Float64}`) adc kspace
"""
get_kspace(hoseq::HO_Sequence; Δt=1,
skip_rf=zeros(Bool, sum(is_RF_on.(hoseq.SEQ)))) = begin
	t, Δt = get_variable_times(hoseq.SEQ; Δt)
	
	Gx, Gy, Gz = get_grads(hoseq.SEQ, t)
	H0, H1, H2, H3, H4, H5, H6, H7, H8 = get_grads(hoseq, t)
	t = t[1:end-1]

	G_nominal = Dict(1=>Gx, 2=>Gy, 3=>Gz)
	G_skope = Dict(1=>H1, 2=>H2, 3=>H3)
	
	K_nominal, K_nominal_adc = grad2kspace(hoseq.SEQ, G_nominal, t, Δt, skip_rf)
	K_skope, K_skope_adc = grad2kspace(hoseq.SEQ, G_skope, t, Δt, skip_rf)
	return K_nominal, K_nominal_adc, K_skope, K_skope_adc
end

function grad2kspace(seq::Sequence, G, t, Δt, skip_rf)
	#kspace
	Nt = length(t)
	k = zeros(Nt,3)
	#get_RF_center_breaks
	idx_rf, rf_type = KomaMRIBase.get_RF_types(seq, t)
	parts = kfoldperm(Nt, 1; breaks=idx_rf)
	for i = 1:3
		kf = 0
		for (rf, p) in enumerate(parts)
			k[p,i] = cumtrapz(Δt[p]', G[i][p[1]:p[end]+1]')[:] #This is the exact integral
			if rf > 1
				if !skip_rf[rf-1]
					if rf_type[rf-1] == 0 # Excitation
						k[p,i] .-= 0
					elseif rf_type[rf-1] == 1 # Refocuse
						k[p,i] .-= kf
					end
				else
					k[p,i] .+= kf
				end
			end
			kf = k[p[end],i]
		end
	end
	kspace = γ * k #[m^-1]
	#Interp, as Gradients are generally piece-wise linear, the integral is piece-wise quadratic
	#Nevertheless, the integral is sampled at the ADC times so a linear interp is sufficient
	#TODO: check if this interpolation is necessary
	ts = t .+ Δt
	t_adc =  get_adc_sampling_times(seq)
	kx_adc = linear_interpolation(ts,kspace[:,1],extrapolation_bc=0)(t_adc)
	ky_adc = linear_interpolation(ts,kspace[:,2],extrapolation_bc=0)(t_adc)
	kz_adc = linear_interpolation(ts,kspace[:,3],extrapolation_bc=0)(t_adc)
	kspace_adc = [kx_adc ky_adc kz_adc]
	#Final
	kspace, kspace_adc	
end 
