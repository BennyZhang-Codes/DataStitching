import KomaMRI.KomaMRICore: run_spin_precession!, run_spin_excitation!
import KomaMRI.KomaMRICore: output_Ndim, initialize_spins_state, Bloch, sim_output_dim

Base.@kwdef struct BlochHighOrder <: SimulationMethod 
    skope::Bool = true
    ho0::Bool = true
    ho1::Bool = true
    ho2::Bool = true
end

export BlochHighOrder
Base.show(io::IO, b::BlochHighOrder) = begin
	print(io, "BlochHighOrder(\n      zero order=$(b.ho0)\n     first order=$(b.ho1)\n    second order=$(b.ho2)\n)")
end

output_Ndim(sim_method::BlochHighOrder) =  output_Ndim(Bloch())#time-points x coils

function sim_output_dim(obj::Phantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochHighOrder) where {T<:Real}
    return sim_output_dim(obj, seq, sys, Bloch()) 
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::Phantom{T}, sim_method::BlochHighOrder) where {T<:Real}
    return initialize_spins_state(obj, Bloch())
end

function run_spin_precession!(p::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochHighOrder) where {T<:Real}
    seq = hoseqd.seqd
    #Simulation
    #Motion
    xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
    #Effective field
    Bzh0 = hoseqd.h0'
    Bzh1 = xt .* hoseqd.h1'
    Bzh2 = yt .* hoseqd.h2'
    Bzh3 = zt .* hoseqd.h3'
    Bzh4 = (xt .* yt) .* hoseqd.h4'
    Bzh5 = (zt .* yt) .* hoseqd.h5'
    Bzh6 = (3zt.^2-(xt.^2 .+ yt.^2 .+ zt.^2)) .* hoseqd.h6'
    Bzh7 = (xt .* zt) .* hoseqd.h7'
    Bzh8 = (xt.^2 .+ yt.^2) .* hoseqd.h8'

    Bzho1 = Bzh1 .+ Bzh2 .+ Bzh3
    Bzho2 = Bzh4 .+ Bzh5 .+ Bzh6 .+ Bzh7 .+ Bzh8

    Bz0 = sim_method.ho0 ? Bzh0 : 0
    Bz1 = sim_method.ho1 ? Bzho1 : xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz'
    Bz2 = sim_method.ho2 ? Bzho2 : 0

    Bz = Bz0 .+ Bz1 .+ Bz2 .+ p.Δw / T(2π * γ)
    #Rotation
    if is_ADC_on(seq)
        ϕ = T(-2π * γ) .* cumtrapz(seq.Δt', Bz)
    else
        ϕ = T(-2π * γ) .* trapz(seq.Δt', Bz)
    end
    #Mxy precession and relaxation, and Mz relaxation
    tp = cumsum(seq.Δt) # t' = t - t0
    dur = sum(seq.Δt)   # Total length, used for signal relaxation
    Mxy = [M.xy M.xy .* exp.(1im .* ϕ .- tp' ./ p.T2)] #This assumes Δw and T2 are constant in time
    M.xy .= Mxy[:, end]
    M.z  .= M.z .* exp.(-dur ./ p.T1) .+ p.ρ .* (1 .- exp.(-dur ./ p.T1))
    #Acquired signal
    sig .= transpose(sum(Mxy[:, findall(seq.ADC)]; dims=1)) #<--- TODO: add coil sensitivities
    return nothing
end





function run_spin_excitation!(p::Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochHighOrder) where {T<:Real}
    return KomaMRICore.run_spin_excitation!(p, hoseqd.seqd, sig, M, Bloch())
end