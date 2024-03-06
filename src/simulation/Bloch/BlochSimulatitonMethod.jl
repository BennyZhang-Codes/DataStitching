import KomaMRI.KomaMRICore: run_spin_precession!, run_spin_excitation!

"""
    run_spin_precession(obj, hoseqd, Xt, sig)

Simulates an MRI sequence `seq` on the Phantom `obj` for time points `t`. It calculates S(t)
= ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). It performs the simulation in free
precession.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct

# Returns
- `S`: (`Vector{ComplexF64}`) raw signal over time
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_spin_precession!(p::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
    return run_spin_precession!(p, hoseqd.seqd, sig, M, sim_method)
end

"""
    M0 = run_spin_excitation(obj, hoseqd, M0)

It gives rise to a rotation of `M0` with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis.

# Arguments
- `obj`: (`::Phantom`) Phantom struct (actually, it's a part of the complete phantom)
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (actually, it's
    a part of the complete Mag vector and it's a part of the initial state for the next
    precession simulation step)
"""
function run_spin_excitation!(p::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
    return run_spin_excitation!(p, hoseqd.seqd, sig, M, sim_method)
end
