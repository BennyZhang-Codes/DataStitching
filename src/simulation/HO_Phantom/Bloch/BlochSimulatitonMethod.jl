

"""
    run_spin_precession(obj, hoseqd, Xt, sig)

Simulates an MRI sequence `seq` on the HO_Phantom `obj` for time points `t`. It calculates S(t)
= ∑ᵢ ρ(xᵢ) exp(- t/T2(xᵢ) ) exp(- 𝒊 γ ∫ Bz(xᵢ,t)). It performs the simulation in free
precession.

# Arguments
- `obj`: (`::HO_Phantom`) HO_Phantom struct (actually, it's a part of the complete phantom)
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct

# Returns
- `S`: (`Vector{ComplexF64}`) raw signal over time
- `M0`: (`::Vector{Mag}`) final state of the Mag vector
"""
function run_spin_precession!(p::HO_Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
    seq = hoseqd.seqd
    #Simulation
    #Motion
    xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
    #Effective field
    Bz = xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz' .+ p.Δw / T(2π * γ)
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

"""
    M0 = run_spin_excitation(obj, hoseqd, M0)

It gives rise to a rotation of `M0` with an angle given by the efective magnetic field
(including B1, gradients and off resonance) and with respect to a rotation axis.

# Arguments
- `obj`: (`::HO_Phantom`) HO_Phantom struct (actually, it's a part of the complete phantom)
- `hoseqd`: (`::HO_DiscreteSequence`) HO_DiscreteSequence struct

# Returns
- `M0`: (`::Vector{Mag}`) final state of the Mag vector after a rotation (actually, it's
    a part of the complete Mag vector and it's a part of the initial state for the next
    precession simulation step)
"""
function run_spin_excitation!(p::HO_Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::Bloch) where {T<:Real}
    seq = hoseqd.seqd
    #Simulation
    for s ∈ seq #This iterates over seq, "s = seq[i,:]"
        #Motion
        xt = p.x .+ p.ux(p.x, p.y, p.z, s.t)
        yt = p.y .+ p.uy(p.x, p.y, p.z, s.t)
        zt = p.z .+ p.uz(p.x, p.y, p.z, s.t)
        #Effective field
        ΔBz = p.Δw ./ T(2π * γ) .- s.Δf ./ T(γ) # ΔB_0 = (B_0 - ω_rf/γ), Need to add a component here to model scanner's dB0(xt,yt,zt)
        Bz = (s.Gx .* xt .+ s.Gy .* yt .+ s.Gz .* zt) .+ ΔBz
        B = sqrt.(abs.(s.B1) .^ 2 .+ abs.(Bz) .^ 2)
        B[B .== 0] .= eps(T)
        #Spinor Rotation
        φ = T(-2π * γ) * (B .* s.Δt) # TODO: Use trapezoidal integration here (?),  this is just Forward Euler
        mul!( Q(φ, s.B1 ./ B, Bz ./ B), M )
        #Relaxation
        M.xy .= M.xy .* exp.(-s.Δt ./ p.T2)
        M.z  .= M.z  .* exp.(-s.Δt ./ p.T1) .+ p.ρ .* (1 .- exp.(-s.Δt ./ p.T1))
    end
    #Acquired signal
    #sig .= -1.4im #<-- This was to test if an ADC point was inside an RF block
    return nothing
end
