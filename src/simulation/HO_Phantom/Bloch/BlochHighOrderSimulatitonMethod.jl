

# output_Ndim(sim_method::BlochHighOrder) =  output_Ndim(Bloch())#time-points x coils

function sim_output_dim(obj::HO_Phantom{T}, seq::Sequence, sys::Scanner, sim_method::BlochHighOrder) where {T<:Real}
    _, nCoil = size(obj.csm)
    return (sum(seq.ADC.N), nCoil) #Nt x Ncoils, This should consider the coil info from sys
end

"""Magnetization initialization for Bloch simulation method."""
function initialize_spins_state(obj::HO_Phantom{T}, sim_method::BlochHighOrder) where {T<:Real}
    Nspins = length(obj)
    Mxy = zeros(T, Nspins)
    Mz = obj.ρ
    Xt = Mag{T}(Mxy, Mz)
    return Xt, obj
end

function run_spin_precession!(p::HO_Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochHighOrder) where {T<:Real}
    seq = hoseqd.seqd
    #Simulation
    #Coil sensitivity
    csm = p.csm  # nSpin x nCoil
    #Motion
    xt = p.x .+ p.ux(p.x, p.y, p.z, seq.t')
    yt = p.y .+ p.uy(p.x, p.y, p.z, seq.t')
    zt = p.z .+ p.uz(p.x, p.y, p.z, seq.t')
    #Effective field
    Bzh1 = xt .* hoseqd.h1'
    Bzh2 = yt .* hoseqd.h2'
    Bzh3 = zt .* hoseqd.h3'
    Bzh4 = (xt .* yt) .* hoseqd.h4'
    Bzh5 = (zt .* yt) .* hoseqd.h5'
    Bzh6 = (3zt.^2-(xt.^2 .+ yt.^2 .+ zt.^2)) .* hoseqd.h6'
    Bzh7 = (xt .* zt) .* hoseqd.h7'
    Bzh8 = (xt.^2 .- yt.^2) .* hoseqd.h8'

    Bz0 = sim_method.ho0 ? hoseqd.h0' : 0
    Bz1 = sim_method.ho1 ? Bzh1 .+ Bzh2 .+ Bzh3 : xt .* seq.Gx' .+ yt .* seq.Gy' .+ zt .* seq.Gz'
    Bz2 = sim_method.ho2 ? Bzh4 .+ Bzh5 .+ Bzh6 .+ Bzh7 .+ Bzh8 : 0

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
    # sig .= transpose(sum(Mxy[:, findall(seq.ADC)] .* csm; dims=1)) #<--- TODO: add coil sensitivities
    sig .= transpose(Mxy[:, findall(seq.ADC)]) * csm  # (nADC x nSpin) * (nSpin x nCoil)
    return nothing
end





function run_spin_excitation!(p::HO_Phantom{T}, hoseqd::HO_DiscreteSequence{T}, sig::AbstractArray{Complex{T}},
    M::Mag{T}, sim_method::BlochHighOrder) where {T<:Real}
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