export SignalOp
# recon
using LinearOperators

mutable struct SignalOp{T,F1,F2} <: AbstractLinearOperator{T}
  nrow :: Int
  ncol :: Int
  symmetric :: Bool
  hermitian :: Bool
  prod! :: Function
  tprod! :: F1
  ctprod! :: F2
  nprod :: Int
  ntprod :: Int
  nctprod :: Int
  args5 :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5 :: Vector{T}
  Mtu5 :: Vector{T}
end

LinearOperators.storage_type(op::SignalOp) = typeof(op.Mv5)

"""
    SignalOp(shape::NTuple{D,Int64}, tr::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D

generates a `SignalOp` which explicitely evaluates the MRI Fourier signal encoding operator.

# Arguments:
* `shape::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr::Trajectory`          - Trajectory with the kspace nodes to sample
                              this results in complex valued image even for real-valued input.
* `Nblocks`                 - split trajectory into `Nblocks` blocks to avoid memory overflow.
* `use_gpu`                 - use GPU for signal encoding/decoding(default: `true`).
"""
function SignalOp(shape::NTuple{D,Int64}, tr::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D
    nodes = Float64.(kspaceNodes(tr))
    times = readoutTimes(tr)
    nrow = size(nodes,2)
    ncol = prod(shape)
    k = size(nodes,2) # number of nodes
    Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k
    
    n = k÷Nblocks # number of nodes per block
    parts = [n for i=1:Nblocks] # number of nodes per block
    parts = [1+n*(i-1):n*i for i=1:Nblocks]
    if k%Nblocks!= 0
        push!(parts, n*Nblocks+1:k)
    end

    Nx, Ny = shape
    x, y = 1:Nx, 1:Ny
    x, y = vec(x .+ y'*0), vec(x*0 .+ y') #grid points
    x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1

    @info "SignalOp Nblocks=$Nblocks, use_gpu=$use_gpu"
    return SignalOp{ComplexF64,Nothing,Function}(nrow, ncol, false, false
                , (res,xm)->(res .= prod_SignalOp(xm, x, y, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu))
                , nothing
                , (res,ym)->(res .= ctprod_SignalOp(ym, x, y, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu))
                , 0,0,0, false, false, false, ComplexF64[], ComplexF64[])
end

function prod_SignalOp(xm::Vector{T}, x::Vector{Float64}, y::Vector{Float64}, nodes::Matrix{Float64};
    Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
    @info "SignalOp prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    out = zeros(ComplexF64,size(nodes,2))
    if use_gpu
        out = out |> gpu
        nodes = nodes |> gpu
        xm = xm |> gpu
        x = x |> gpu
        y = y |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        kx, ky = nodes[1,p], nodes[2,p]
        ϕ = (kx .* x') .+ (ky .* y')
        e = exp.(-2*1im*pi*ϕ)
        out[p] =  e * xm
        next!(progress_bar, showvalues=[(:Nblocks, block)])
    end
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function ctprod_SignalOp(xm::Vector{T}, x::Vector{Float64}, y::Vector{Float64}, nodes::Matrix{Float64}; 
    Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
    @info "SignalOp ctprod Nblocks=$Nblocks, use_gpu=$use_gpu"
    out = zeros(ComplexF64, size(x, 1))
    if use_gpu
        out = out |> gpu
        nodes = nodes |> gpu
        xm = xm |> gpu
        x = x |> gpu
        y = y |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        kx, ky = nodes[1,p], nodes[2,p]
        ϕ = (x .* kx') .+ (y .* ky')
        e = exp.(-2*1im*pi*ϕ)
        out +=  conj(e) * xm[p]
        next!(progress_bar, showvalues=[(:Nblocks, block)])
    end
    if use_gpu
        out = out |> cpu
    end
  return vec(out)
end

function prod_SignalOp(xm::SubArray{T}, x::Vector{Float64}, y::Vector{Float64}, nodes::Matrix{Float64};
    Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
	xm = xm[:,1]
  return prod_SignalOp(xm, x, y, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu)
end

function ctprod_SignalOp(xm::SubArray{T}, x::Vector{Float64}, y::Vector{Float64}, nodes::Matrix{Float64}; 
    Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
	xm = xm[:,1]
  return ctprod_SignalOp(xm, x, y, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu)
end


function Base.adjoint(op::SignalOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
