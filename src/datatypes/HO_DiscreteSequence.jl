import KomaMRI.KomaMRIBase: discretize, get_rfs, default_sampling_params, get_variable_times

"""
    seqd = HO_DiscreteSequence(seqd, h0, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15)

# Description
    A sampled version of a HO_Sequence struct, containing vectors for event amplitudes at specified
times. DiscreteSequence is the struct used for simulation.

# Arguments
- `seqd`: (`::DiscreteSequence`) 
- `h0`: (`::AbstractVector{T<:Real}`, `[T]`)    `1`
- `h1`: (`::AbstractVector{T<:Real}`, `[T/m]`)  `x`
- `h2`: (`::AbstractVector{T<:Real}`, `[T/m]`)  `y`
- `h3`: (`::AbstractVector{T<:Real}`, `[T/m]`)  `z`
- `h4`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `xy`
- `h5`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `zy`
- `h6`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `3z² - (x²	+ y² + z²)`
- `h7`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `xz`
- `h8`: (`::AbstractVector{T<:Real}`, `[T/m²]`) `x² - y²`
- `h9`: (`::AbstractVector{T<:Real}`, `[T/m³]`)  `3yx² - y³`
- `h10`: (`::AbstractVector{T<:Real}`, `[T/m³]`) `xzy`
- `h11`: (`::AbstractVector{T<:Real}`, `[T/m³]`) `5z² - (x² + y² + z²))y`
- `h12`: (`::AbstractVector{T<:Real}`, `[T/m³]`) `5z³ - 3z(x² + y² + z²)`
- `h13`: (`::AbstractVector{T<:Real}`, `[T/m³]`) `5z² - (x² + y² + z²))x`
- `h14`: (`::AbstractVector{T<:Real}`, `[T/m³]`) `x²z - y²z`
- `h15`: (`::AbstractVector{T<:Real}`, `[T/m³]`) `x³ - 3xy²`
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
    h9::AbstractVector{T}
    h10::AbstractVector{T}
    h11::AbstractVector{T}
    h12::AbstractVector{T}
    h13::AbstractVector{T}
    h14::AbstractVector{T}
    h15::AbstractVector{T}
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
                        hoseqd.h8[i, :],
                        hoseqd.h9[i, :],
                        hoseqd.h10[i, :],
                        hoseqd.h11[i, :],
                        hoseqd.h12[i, :],
                        hoseqd.h13[i, :],
                        hoseqd.h14[i, :],
                        hoseqd.h15[i, :]
                        )
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
                        hoseqd.h8[i.start:i.stop+1],
                        hoseqd.h9[i.start:i.stop+1],
                        hoseqd.h10[i.start:i.stop+1],
                        hoseqd.h11[i.start:i.stop+1],
                        hoseqd.h12[i.start:i.stop+1],
                        hoseqd.h13[i.start:i.stop+1],
                        hoseqd.h14[i.start:i.stop+1],
                        hoseqd.h15[i.start:i.stop+1]
                        )
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
                                hoseqd.h8[i.start:i.stop+1],
                                hoseqd.h9[i.start:i.stop+1],
                                hoseqd.h10[i.start:i.stop+1],
                                hoseqd.h11[i.start:i.stop+1],
                                hoseqd.h12[i.start:i.stop+1],
                                hoseqd.h13[i.start:i.stop+1],
                                hoseqd.h14[i.start:i.stop+1],
                                hoseqd.h15[i.start:i.stop+1]
                                )
end
Base.iterate(hoseqd::HO_DiscreteSequence) = (hoseqd[1], 2)
Base.iterate(hoseqd::HO_DiscreteSequence, i) = (i <= length(hoseqd)) ? (hoseqd[i], i+1) : nothing

"""
    hoseqd = discretize(hoseq::HO_Sequence; sampling_params=default_sampling_params())

# Description
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
function discretize(hoseq::HO_Sequence; sampling_params=default_sampling_params())::HO_DiscreteSequence
    t, Δt      = get_variable_times(hoseq.SEQ; Δt=sampling_params["Δt"], Δt_rf=sampling_params["Δt_rf"])
    B1, Δf     = get_rfs(hoseq.SEQ, t)
    Gx, Gy, Gz = get_grads(hoseq.SEQ, t)
    tadc       = get_adc_sampling_times(hoseq.SEQ)
    ADCflag    = [any(tt .== tadc) for tt in t] #Displaced 1 dt, sig[i]=S(ti+dt)
    seqd       = DiscreteSequence(Gx, Gy, Gz, complex.(B1), Δf, ADCflag, t, Δt)
    H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15 = get_grads(hoseq, t)
    hoseqd = HO_DiscreteSequence(seqd, H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15)
    return hoseqd
end



