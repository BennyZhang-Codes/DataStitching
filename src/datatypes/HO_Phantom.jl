"""
    obj = HO_Phantom(name, x, y, z, ρ, T1, T2, T2s, Δw, Dλ1, Dλ2, Dθ, ux, uy, uz, csm)

The HO_Phantom struct. Most of its field names are vectors, with each element associated with
a property value representing a spin. This struct serves as an input for the simulation.

# Arguments
- `name`: (`::String`) phantom name
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `ρ`: (`::AbstractVector{T<:Real}`) spin proton density vector
- `T1`: (`::AbstractVector{T<:Real}`, `[s]`) spin T1 parameter vector
- `T2`: (`::AbstractVector{T<:Real}`, `[s]`) spin T2 parameter vector
- `T2s`: (`::AbstractVector{T<:Real}`, `[s]`) spin T2s parameter vector
- `Δw`: (`::AbstractVector{T<:Real}`, `[rad/s]`) spin off-resonance parameter vector
- `Dλ1`: (`::AbstractVector{T<:Real}`) spin Dλ1 (diffusion) parameter vector
- `Dλ2`: (`::AbstractVector{T<:Real}`) spin Dλ2 (diffusion) parameter vector
- `Dθ`: (`::AbstractVector{T<:Real}`) spin Dθ (diffusion) parameter vector
- `ux`: (`::Function`) displacement field in the x-axis
- `uy`: (`::Function`) displacement field in the y-axis
- `uz`: (`::Function`) displacement field in the z-axis
- `nCoil`: (`::Int64`) number of coils in the phantom
- `csm`: (`::AbstractArray{Complex{T}, 2}`) Coil-Sensitivity Map (CSM) matrix

# Returns
- `obj`: (`::HO_Phantom`) HO_Phantom struct

# Examples
```julia-repl
julia> obj = HO_Phantom(x=[0.0])

julia> obj.ρ
```
"""
 @with_kw mutable struct HO_Phantom{T<:Real}
    name::String = "spins"
	x::AbstractVector{T}
	y::AbstractVector{T} = zeros(size(x))
	z::AbstractVector{T} = zeros(size(x))
	ρ::AbstractVector{T} = ones(size(x))
	T1::AbstractVector{T} = ones(size(x)) * 1_000_000
	T2::AbstractVector{T} = ones(size(x)) * 1_000_000
	T2s::AbstractVector{T} = ones(size(x)) * 1_000_000
	#Off-resonance related
	Δw::AbstractVector{T} = zeros(size(x))
	#χ::Vector{SusceptibilityModel}
	#Diffusion
	Dλ1::AbstractVector{T} = zeros(size(x))
	Dλ2::AbstractVector{T} = zeros(size(x))
	Dθ::AbstractVector{T} =  zeros(size(x))
	#Diff::Vector{DiffusionModel}  #Diffusion map
	#Motion
	ux::Function = (x,y,z,t)->0
	uy::Function = (x,y,z,t)->0
	uz::Function = (x,y,z,t)->0
    #Coil-Sensitivity
    csm::AbstractArray{Complex{T}, 2} = zeros(size(x), 1)
end

"""Size and length of a phantom"""
size(x::HO_Phantom) = size(x.ρ)
Base.length(x::HO_Phantom) = length(x.ρ)
# To enable to iterate and broadcast over the HO_Phantom
Base.iterate(x::HO_Phantom) = (x[1], 2)
Base.iterate(x::HO_Phantom, i::Integer) = (i <= length(x)) ? (x[i], i+1) : nothing
Base.lastindex(x::HO_Phantom) = length(x)
Base.getindex(x::HO_Phantom, i::Integer) = x[i:i]

"""Compare two phantoms"""
Base.isapprox(obj1::HO_Phantom, obj2::HO_Phantom)  = begin
    obj1.x     ≈ obj2.x    &&
    obj1.y     ≈ obj2.y    &&
    obj1.z     ≈ obj2.z    &&
    obj1.ρ     ≈ obj2.ρ    &&
    obj1.T1    ≈ obj2.T1   &&
    obj1.T2    ≈ obj2.T2   &&
    obj1.T2s   ≈ obj2.T2s  &&
    obj1.Δw    ≈ obj2.Δw   &&
    obj1.Dλ1   ≈ obj2.Dλ1  &&
    obj1.Dλ2   ≈ obj2.Dλ2  &&
    obj1.Dθ    ≈ obj2.Dθ   &&
	obj1.csm  ≈ obj2.csm 
end

"""
Separate object spins in a sub-group
"""
Base.getindex(obj::HO_Phantom, p::AbstractRange) = begin
	HO_Phantom(name=obj.name,
			x=obj.x[p],
			y=obj.y[p],
			z=obj.z[p],
			ρ=obj.ρ[p],
			T1=obj.T1[p],
			T2=obj.T2[p],
			T2s=obj.T2s[p],
			Δw=obj.Δw[p],
			#Diff=obj.Diff[p], #TODO!
			Dλ1=obj.Dλ1[p],
			Dλ2=obj.Dλ2[p],
			Dθ=obj.Dθ[p],
			#Χ=obj.Χ[p], #TODO!
			ux=obj.ux,
			uy=obj.uy,
			uz=obj.uz,
			csm=obj.csm[p,:]
			)
end

"""Separate object spins in a sub-group (lightweigth)."""
Base.view(obj::HO_Phantom, p::AbstractRange) = begin
	@views HO_Phantom(name=obj.name,
			x=obj.x[p],
			y=obj.y[p],
			z=obj.z[p],
			ρ=obj.ρ[p],
			T1=obj.T1[p],
			T2=obj.T2[p],
			T2s=obj.T2s[p],
			Δw=obj.Δw[p],
			#Diff=obj.Diff[p], #TODO!
			Dλ1=obj.Dλ1[p],
			Dλ2=obj.Dλ2[p],
			Dθ=obj.Dθ[p],
			#Χ=obj.Χ[p], #TODO!
			ux=obj.ux,
			uy=obj.uy,
			uz=obj.uz,
			csm=obj.csm[p,:]
			)
end

"""Addition of phantoms"""
+(s1::HO_Phantom,s2::HO_Phantom) = begin
	HO_Phantom(name=s1.name*"+"*s2.name,
		x=[s1.x;s2.x],
		y=[s1.y;s2.y],
		z=[s1.z;s2.z],
		ρ=[s1.ρ;s2.ρ],
		T1=[s1.T1;s2.T1],
		T2=[s1.T2;s2.T2],
		T2s=[s1.T2s;s2.T2s],
		Δw=[s1.Δw;s2.Δw],
		#Diff=obj.Diff[p], #TODO!
		Dλ1=[s1.Dλ1;s2.Dλ1],
		Dλ2=[s1.Dλ2;s2.Dλ2],
		Dθ=[s1.Dθ;s2.Dθ],
		#Χ=obj.Χ[p], #TODO!
		ux=s1.ux,
		uy=s1.uy,
		uz=s1.uz,
		csm=[s1.csm;s2.csm]
	)
end

"""Scalar multiplication of a phantom"""
*(α::Real,obj::HO_Phantom) = begin
	HO_Phantom(name=obj.name,
		x=obj.x,
		y=obj.y,
		z=obj.z,
		ρ=α*obj.ρ, #Only affects the proton density
		T1=obj.T1,
		T2=obj.T2,
		T2s=obj.T2s,
		Δw=obj.Δw,
		#Diff=obj.Diff[p], #TODO!
		Dλ1=obj.Dλ1,
		Dλ2=obj.Dλ2,
		Dθ=obj.Dθ,
		#Χ=obj.Χ[p], #TODO!
		ux=obj.ux,
		uy=obj.uy,
		uz=obj.uz,
		csm=obj.csm
	)
end

