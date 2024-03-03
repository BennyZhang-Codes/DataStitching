"""
    seqd = HO_DiscreteSequence(Gx, Gy, Gz, B1, Δf, ADC)

A sampled version of a HO_Sequence struct, containing vectors for event amplitudes at specified
times. DiscreteSequence is the struct used for simulation.

# Arguments
- `seqd`: (`::DiscreteSequence`) 
- `h0`: (`::AbstractVector{T<:Real}`, `[T]`) h0-gradient vector
- `h1`: (`::AbstractVector{T<:Real}`, `[T/m]`) h1-gradient vector
- `h2`: (`::AbstractVector{T<:Real}`, `[T/m]`) h2-gradient vector
- `h3`: (`::AbstractVector{T<:Real}`, `[T/m]`) h3-gradient vector
- `h4`: (`::AbstractVector{T<:Real}`, `[T/m²]`) h4-gradient vector
- `h5`: (`::AbstractVector{T<:Real}`, `[T/m²]`) h5-gradient vector
- `h6`: (`::AbstractVector{T<:Real}`, `[T/m²]`) h6-gradient vector
- `h7`: (`::AbstractVector{T<:Real}`, `[T/m²]`) h7-gradient vector
- `h8`: (`::AbstractVector{T<:Real}`, `[T/m²]`) h8-gradient vector

# Returns
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct
"""
struct HO_DiscreteSequence{T<:Real}
    seqd::DiscreteSequence
    h0::AbstractVector{T}
    h1::AbstractVector{T}
    h2::AbstractVector{T}
    h3::AbstractVector{T}
    h4::AbstractVector{T}
    h5::AbstractVector{T}
    h6::AbstractVector{T}
    h7::AbstractVector{T}
    h8::AbstractVector{T}
end

"""
    hoseqd = HO_discretize(hoseq::HO_Sequence; sampling_params=default_sampling_params())

This function returns a sampled Sequence struct with RF and gradient time refinements
based on simulation parameters.

# Arguments
- `hoseq`: (`::HO_Sequence`) sequence

# Keywords
- `sampling_params`: (`::Dict{String, Any}`, `=default_sampling_params()`) sampling
    parameter dictionary

# Returns
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct
"""
function HO_discretize(hoseq::HO_Sequence; sampling_params=KomaMRIBase.default_sampling_params())
    t, Δt      = KomaMRIBase.get_variable_times(hoseq.SEQ; Δt=sampling_params["Δt"], Δt_rf=sampling_params["Δt_rf"])
    B1, Δf     = KomaMRIBase.get_rfs(hoseq.SEQ, t)
    Gx, Gy, Gz = KomaMRIBase.get_grads(hoseq.SEQ, t)
    tadc       = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ)
    ADCflag    = [any(tt .== tadc) for tt in t] #Displaced 1 dt, sig[i]=S(ti+dt)
    seqd       = KomaMRIBase.DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    H0, H1, H2, H3, H4, H5, H6, H7, H8 = HO_get_grads(hoseq, t)
    hoseqd = HO_DiscreteSequence(seqd, H0, H1, H2, H3, H4, H5, H6, H7, H8)
    return hoseqd
end



function HO_get_theo_Gi(hoseq::HO_Sequence, idx)
    N = length(hoseq.SEQ)
    T0 = KomaMRIBase.get_block_start_times(hoseq.SEQ)
    t = vcat([KomaMRIBase.get_theo_t(hoseq.GR_skope[idx,i]) .+ T0[i] for i=1:N]...)
    G = vcat([KomaMRIBase.get_theo_A(hoseq.GR_skope[idx,i]) for i=1:N]...) #; off_val=0 <---potential solution
	Interpolations.deduplicate_knots!(t; move_knots=true)
	return (t, G)
end


function HO_get_grads(hoseq::HO_Sequence, t::Vector)
    h0 = HO_get_theo_Gi(hoseq, 1)
    h1 = HO_get_theo_Gi(hoseq, 2)
    h2 = HO_get_theo_Gi(hoseq, 3)
    h3 = HO_get_theo_Gi(hoseq, 4)
    h4 = HO_get_theo_Gi(hoseq, 5)
    h5 = HO_get_theo_Gi(hoseq, 6)
    h6 = HO_get_theo_Gi(hoseq, 7)
    h7 = HO_get_theo_Gi(hoseq, 8)
    h8 = HO_get_theo_Gi(hoseq, 9)
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