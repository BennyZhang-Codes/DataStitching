struct Bloch_HO1 <: SimulationMethod end

export Bloch_HO1


function run_spin_precession!(p::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch_HO1) where {T<:Real}
    return KomaMRICore.run_spin_precession!(p, hoseqd.seqd, sig, M, sim_method)
end





function run_spin_excitation!(p::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch_HO1) where {T<:Real}
    return KomaMRICore.run_spin_excitation!(p, hoseqd.seqd, sig, M, sim_method)
end