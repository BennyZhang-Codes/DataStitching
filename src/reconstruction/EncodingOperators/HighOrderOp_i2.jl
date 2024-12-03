export HighOrderOp_i2
# recon
using LinearOperators

mutable struct HighOrderOp_i2{T,F1,F2} <: AbstractLinearOperator{T}
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

LinearOperators.storage_type(op::HighOrderOp_i2) = typeof(op.Mv5)

"""
    HighOrderOp_i2(reconsize::NTuple{D,Int64}, tr_nominal::Trajectory, tr_measured::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D

# Description
    generates a `HighOrderOp_i2` which explicitely evaluates the MRI Fourier HighOrder encoding operator.

# Arguments:
* `reconsize::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr_nominal::Trajectory`  - (must be 3 rows [x, y, z]') nominal Trajectory without normalization.
* `tr_measured::Trajectory`    - (must be 9 rows [h0, h1, h2, h3, h4, h5, h6, h7, h8]') measured Trajectory without normalization. 
``
                              this results in complex valued image even for real-valued input.
* `Nblocks`                 - split trajectory into `Nblocks` blocks to avoid memory overflow.
* `use_gpu`                 - use GPU for HighOrder encoding/decoding(default: `true`).
"""
function HighOrderOp_i2(
    reconsize::NTuple{D,Int64}, 
    tr_nominal::Trajectory, 
    tr_measured::Trajectory, 
    sim_method::BlochHighOrder; 
    fieldmap::Matrix{T}=zeros(Float64,reconsize), 
    csm::Array{Complex{T}, 3}=ones(ComplexF64,reconsize..., 1), 
    Nblocks::Int64=50, 
    use_gpu::Bool=true, 
    verbose::Bool=false, 
    Δx::Float64=1e-3, 
    Δy::Float64=1e-3,
    grid::Int64=1,
    ) where {D, T<:AbstractFloat}

    nodes_measured = T.(kspaceNodes(tr_measured))
    nodes_nominal = T.(kspaceNodes(tr_nominal))
    times = T.(readoutTimes(tr_measured))
    @assert size(nodes_measured,1) == 9 "nodes for measured must have 9 rows"
    @assert size(nodes_nominal,1) == 3 "nodes for nominal must have 3 rows"
    @assert size(fieldmap) == reconsize "fieldmap must have same size as reconsize"
    @assert size(csm)[1:2] == reconsize "fieldmap must have same size as reconsize"

    nX, nY = reconsize
    nCha = size(csm, 3)
    nSam = size(nodes_measured,2)
    nRow = nSam * nCha
    nCol = nVox = prod(reconsize)
    
    csm = reshape(csm, nX*nY, nCha)  # [nX*nY, nCha]
    fieldmap = vec(fieldmap)

    k = nSam # number of nodes
    Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k
    
    n = k÷Nblocks # number of nodes per block
    parts = [n for i=1:Nblocks] # number of nodes per block
    parts = [1+n*(i-1):n*i for i=1:Nblocks]
    if k%Nblocks!= 0
        push!(parts, n*Nblocks+1:k)
    end

    
    x, y = 1:nX, 1:nY
    if grid == 1      # x up->down, y left->right
        x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- nX/2 .- 1, y .- nY/2 .- 1
    elseif grid == 2
        x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (nX+1)/2, y .- (nY+1)/2
    elseif grid == 3  # x left->right, y down->up
        x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (nX+1)/2, y .- (nY+1)/2
    elseif grid == 4  # x left->right, y up->down
        x, y, z = vec(x*0.0 .+ y'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (nX+1)/2, y .- (nY+1)/2
    elseif grid == 5  # x right->left, y up->down
        x, y, z = vec(x*0.0 .+ y[end:-1:1]'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (nX+1)/2, y .- (nY+1)/2
    end
    # print(x)
    x, y = x * Δx, y * Δy 


    @info "HighOrderOp_i2 Nblocks=$Nblocks, use_gpu=$use_gpu, Δx=$Δx, Δy=$Δy"
    @info sim_method
    return HighOrderOp_i2{ComplexF64,Nothing,Function}(nRow, nCol, false, false
                , (res,xm)->(res .= prod_HighOrderOp_i2(xm, x, y, z, nVox, nSam, nCha, nodes_measured, nodes_nominal, times, fieldmap, csm;
                                        sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
                , nothing
                , (res,ym)->(res .= ctprod_HighOrderOp_i2(ym, x, y, z, nVox, nSam, nCha, nodes_measured, nodes_nominal, times, fieldmap, csm;
                                        sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
                , 0,0,0, false, false, false, ComplexF64[], ComplexF64[])
end

function prod_HighOrderOp_i2(
    xm::AbstractVector{T}, 
    x::Vector{Float64}, 
    y::Vector{Float64}, 
    z::Vector{Float64}, 
    nVox::Int64, 
    nSam::Int64, 
    nCha::Int64,
    nodes_measured::Matrix{Float64}, 
    nodes_nominal::Matrix{Float64},
    times::Vector{Float64}, 
    fieldmap::Vector{Float64},
    csm::Array{ComplexF64, 2};
    sim_method::BlochHighOrder=BlochHighOrder("111"), 
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where T<:Union{Real,Complex}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp_i2 prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nSam, nCha)
    x0 = ones(Float64, nVox)
    if use_gpu
        out = out |> gpu
        nodes_measured = nodes_measured |> gpu
        nodes_nominal = nodes_nominal |> gpu
        xm = xm |> gpu
        x0 = x0 |> gpu
        x = x |> gpu
        y = y |> gpu
        z = z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csm = csm |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                       @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
        hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
        ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
        ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
        ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
                h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
        ϕB0 = times[p] .* fieldmap'
        ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
        e = exp.(-2*1im*pi*ϕ)
        out[p, :] =  e * (xm .* csm)
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function ctprod_HighOrderOp_i2(
    ym::AbstractVector{T}, 
    x::Vector{Float64}, 
    y::Vector{Float64}, 
    z::Vector{Float64},
    nVox::Int64, 
    nSam::Int64, 
    nCha::Int64,
    nodes_measured::Matrix{Float64}, 
    nodes_nominal::Matrix{Float64},
    times::Vector{Float64}, 
    fieldmap::Vector{Float64},
    csm::Array{ComplexF64, 2};
    sim_method::BlochHighOrder=BlochHighOrder("111"), 
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where T<:Union{Real,Complex}
    csmC = conj.(csm)
    ym = reshape(ym, nSam, nCha)
    # ym = Vector(ym[1:nSam])

    if verbose
        @info "HighOrderOp_i2 ctprod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nVox, nCha)
    x0 = ones(Float64, nVox)

    if use_gpu
        out = out |> gpu
        nodes_measured = nodes_measured |> gpu
        nodes_nominal = nodes_nominal |> gpu
        ym = ym |> gpu
        x0 = x0 |> gpu
        x = x |> gpu
        y = y |> gpu
        z = z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csmC = csmC |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                       @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
        hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
        ϕ0 = sim_method.ho0 ? x0 .* h0' : 0
        ϕ1 = sim_method.ho1 ? (x .* h1') .+ (y .* h2') .+ (z .* h3') : (x .* hx') .+ (y .* hy') .+ (z .* hz')
        ϕ2 = sim_method.ho2 ? (x .* y) .* h4' .+ (z .* y) .* h5' .+ (3z.^2-(x.^2 .+ y.^2 .+ z.^2)) .* h6' .+
            (x .* z) .* h7' .+ (x.^2 .- y.^2) .* h8' : 0
        ϕB0 = fieldmap .* times[p]'
        ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0

        e = exp.(-2*1im*pi*ϕ)
        out +=  conj(e) * ym[p, :]
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
    out = out .* csmC
    if use_gpu
        out = out |> cpu
    end
  return vec(sum(out, dims=2))
end

function Base.adjoint(op::HighOrderOp_i2{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
