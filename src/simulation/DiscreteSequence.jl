import KomaMRI.KomaMRIBase: discretize, get_theo_Gi, get_grads

"""
    seqd = HO_DiscreteSequence(seqd, h0, h1, h2, h3, h4, h5, h6, h7, h8)

A sampled version of a HO_Sequence struct, containing vectors for event amplitudes at specified
times. DiscreteSequence is the struct used for simulation.

# Arguments
- `seqd`: (`::DiscreteSequence`) 
- `h0`: (`::AbstractVector{T<:Real}`, `[T]`)    `1`
- `h1`: (`::AbstractVector{T<:Real}`, `[T/m]`)  `x`
- `h2`: (`::AbstractVector{T<:Real}`, `[T/m]`)  `y`
- `h3`: (`::AbstractVector{T<:Real}`, `[T/m]`)  `z`
- `h4`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `xy`
- `h5`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `yz`
- `h6`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `3z² - (x²	+ y² + z²)`
- `h7`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `xz`
- `h8`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `x² - y²`
# Returns
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct
"""
struct HO_DiscreteSequence{T<:Real}
    seqd::DiscreteSequence{T}
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
@functor HO_DiscreteSequence

Base.length(hoseqd::HO_DiscreteSequence) = length(hoseqd.seqd.Δt)
Base.getindex(hoseqd::HO_DiscreteSequence, i::Integer) = begin
    HO_DiscreteSequence(hoseqd.seqd[i],
                        hoseqd.h0[i, :],
                        hoseqd.h1[i, :],
                        hoseqd.h2[i, :],
                        hoseqd.h3[i, :],
                        hoseqd.h4[i, :],
                        hoseqd.h5[i, :],
                        hoseqd.h6[i, :],
                        hoseqd.h7[i, :],
                        hoseqd.h8[i, :])
end
Base.getindex(hoseqd::HO_DiscreteSequence, i::UnitRange) = begin
    HO_DiscreteSequence(hoseqd.seqd[i.start:i.stop],
                        hoseqd.h0[i.start:i.stop+1],
                        hoseqd.h1[i.start:i.stop+1],
                        hoseqd.h2[i.start:i.stop+1],
                        hoseqd.h3[i.start:i.stop+1],
                        hoseqd.h4[i.start:i.stop+1],
                        hoseqd.h5[i.start:i.stop+1],
                        hoseqd.h6[i.start:i.stop+1],
                        hoseqd.h7[i.start:i.stop+1],
                        hoseqd.h8[i.start:i.stop+1])
end
Base.view(hoseqd::HO_DiscreteSequence, i::UnitRange) = begin
    @views HO_DiscreteSequence(hoseqd.seqd[i.start:i.stop],
                                hoseqd.h0[i.start:i.stop+1],
                                hoseqd.h1[i.start:i.stop+1],
                                hoseqd.h2[i.start:i.stop+1],
                                hoseqd.h3[i.start:i.stop+1],
                                hoseqd.h4[i.start:i.stop+1],
                                hoseqd.h5[i.start:i.stop+1],
                                hoseqd.h6[i.start:i.stop+1],
                                hoseqd.h7[i.start:i.stop+1],
                                hoseqd.h8[i.start:i.stop+1])
end
Base.iterate(hoseqd::HO_DiscreteSequence) = (hoseqd[1], 2)
Base.iterate(hoseqd::HO_DiscreteSequence, i) = (i <= length(hoseqd)) ? (hoseqd[i], i+1) : nothing

"""
    hoseqd = discretize(hoseq::HO_Sequence; sampling_params=default_sampling_params())

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
function discretize(hoseq::HO_Sequence; sampling_params=KomaMRIBase.default_sampling_params())::HO_DiscreteSequence
    t, Δt      = KomaMRIBase.get_variable_times(hoseq.SEQ; Δt=sampling_params["Δt"], Δt_rf=sampling_params["Δt_rf"])
    B1, Δf     = KomaMRIBase.get_rfs(hoseq.SEQ, t)
    Gx, Gy, Gz = KomaMRIBase.get_grads(hoseq.SEQ, t)
    tadc       = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ)
    ADCflag    = [any(tt .== tadc) for tt in t] #Displaced 1 dt, sig[i]=S(ti+dt)
    seqd       = KomaMRIBase.DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    H0, H1, H2, H3, H4, H5, H6, H7, H8 = get_grads(hoseq, t)
    hoseqd = HO_DiscreteSequence(seqd, H0, H1, H2, H3, H4, H5, H6, H7, H8)
    return hoseqd
end



function get_theo_Gi(hoseq::HO_Sequence, idx)
    N = length(hoseq.SEQ)
    T0 = KomaMRIBase.get_block_start_times(hoseq.SEQ)
    t = vcat([KomaMRIBase.get_theo_t(hoseq.GR_skope[idx,i]) .+ T0[i] for i=1:N]...)
    G = vcat([KomaMRIBase.get_theo_A(hoseq.GR_skope[idx,i]) for i=1:N]...) #; off_val=0 <---potential solution
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