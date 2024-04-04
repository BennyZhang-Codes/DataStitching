export HighOrderOp
# recon
using LinearOperators

mutable struct HighOrderOp{T,F1,F2} <: AbstractLinearOperator{T}
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

LinearOperators.storage_type(op::HighOrderOp) = typeof(op.Mv5)

"""
    HighOrderOp(shape::NTuple{D,Int64}, tr::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D

generates a `HighOrderOp` which explicitely evaluates the MRI Fourier HighOrder encoding operator.

# Arguments:
* `shape::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr::Trajectory`          - Trajectory with the kspace nodes to sample
                              this results in complex valued image even for real-valued input.
* `Nblocks`                 - split trajectory into `Nblocks` blocks to avoid memory overflow.
* `use_gpu`                 - use GPU for HighOrder encoding/decoding(default: `true`).
"""
function HighOrderOp(shape::NTuple{D,Int64}, tr::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D
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
    x, y, z = vec(x .+ y'*0), vec(x*0 .+ y'), vec(x*0 .+ y'*0) #grid points
    x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1

    @info "HighOrderOp Nblocks=$Nblocks, use_gpu=$use_gpu"
    return HighOrderOp{ComplexF64,Nothing,Function}(nrow, ncol, false, false
                , (res,xm)->(res .= prod_HighOrderOp(xm, x, y, z, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu))
                , nothing
                , (res,ym)->(res .= ctprod_HighOrderOp(ym, x, y, z, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu))
                , 0,0,0, false, false, false, ComplexF64[], ComplexF64[])
end

function prod_HighOrderOp(xm::Vector{T}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, nodes::Matrix{Float64};
    Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
    @info "HighOrderOp prod Nblocks=$Nblocks, use_gpu=$use_gpu"
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
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = nodes[1,p], nodes[2,p], nodes[3,p], nodes[4,p], nodes[5,p], nodes[6,p], nodes[7,p], nodes[8,p]
        ϕ0 = h0 .* ones(Float64, size(x))
        ϕ1 = (h1 .* x') .+ (h2 .* y') .+ (h3 .* z')
        ϕ2 = h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2)) .+
                h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)'
        ϕ = ϕ0 .+ ϕ1 .+ ϕ2

        # kx, ky = nodes[1,p], nodes[2,p]
        # ϕ = (kx .* x') .+ (ky .* y')

        e = exp.(-2*1im*pi*ϕ)
        out[p] =  e * xm
        next!(progress_bar, showvalues=[(:Nblocks, block)])
    end
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function ctprod_HighOrderOp(xm::Vector{T}, x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, nodes::Matrix{Float64}; 
    Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
    @info "HighOrderOp ctprod Nblocks=$Nblocks, use_gpu=$use_gpu"
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

function Base.adjoint(op::HighOrderOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
