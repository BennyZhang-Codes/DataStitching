export HighOrderOpv2w
# recon
using LinearOperators

mutable struct HighOrderOpv2w{T,F1,F2} <: AbstractLinearOperator{T}
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

LinearOperators.storage_type(op::HighOrderOpv2w) = typeof(op.Mv5)

"""
    HighOrderOpv2w(reconsize::NTuple{D,Int64}, tr_nominal::Trajectory, tr_measured::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D

# Description
    generates a `HighOrderOpv2w` which explicitely evaluates the MRI Fourier HighOrder encoding operator.

# Arguments:
* `reconsize::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr_nominal::Trajectory`  - (must be 3 rows [x, y, z]') nominal Trajectory without normalization.
* `tr_measured::Trajectory`    - (must be 9 rows [h0, h1, h2, h3, h4, h5, h6, h7, h8]') measured Trajectory without normalization. 
``
                              this results in complex valued image even for real-valued input.
* `Nblocks`                 - split trajectory into `Nblocks` blocks to avoid memory overflow.
* `use_gpu`                 - use GPU for HighOrder encoding/decoding(default: `true`).
"""
function HighOrderOpv2w(
    grid::Grid{T},
    tr_kspha::AbstractArray{T, 2}, 
    times::AbstractVector{T}, 
    sim_method::BlochHighOrder; 
    fieldmap::Matrix{T}=zeros(T,(grid.nX, grid.nY)), 
    csm::Array{Complex{T}, 3}=ones(Complex{T},(grid.nX, grid.nY)..., 1), 
    tr_nominal::AbstractArray{T, 2}=tr_kspha[2:4, :],
    Nblocks::Int64=50, 
    use_gpu::Bool=true, 
    verbose::Bool=false, 
    ) where {T<:AbstractFloat}

    @assert size(tr_kspha,  1) == 9 "tr_kspha must have 9 rows for 0th-2nd order terms"
    @assert size(tr_nominal,1) == 3 "tr_nominal must have 3 rows for kx, ky, kz"
    @assert size(fieldmap) == (grid.nX, grid.nY) "FieldMap must have same size as $((grid.nX, grid.nY)) in grid"
    @assert size(csm)[1:2] == (grid.nX, grid.nY) "Coil-SensitivityMap must have same size as $((grid.nX, grid.nY)) in grid"

    nX, nY = grid.nX, grid.nY
    nCha = size(csm, 3)
    nSam = size(tr_kspha,2)
    nRow = nSam * nCha
    nCol = nVox = prod(grid.matrixSize)
    
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

    @info "HighOrderOpv2w Nblocks=$Nblocks, use_gpu=$use_gpu"
    @info grid
    @info sim_method
    return HighOrderOpv2w{Complex{T},Nothing,Function}(nRow, nCol, false, false
                , (res,xm)->(res .= prod_HighOrderOpv2w(xm, grid, nVox, nSam, nCha, tr_kspha, tr_nominal, times, fieldmap, csm;
                                        sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
                , nothing
                , (res,ym)->(res .= ctprod_HighOrderOpv2w(ym, grid, nVox, nSam, nCha, tr_kspha, tr_nominal, times, fieldmap, csm;
                                        sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
                , 0,0,0, false, false, false, Complex{T}[], Complex{T}[])
end

function prod_HighOrderOpv2w(
    xm::AbstractVector{T}, 
    grid::Grid{D},
    nVox::Int64, 
    nSam::Int64, 
    nCha::Int64,
    tr_kspha::AbstractArray{D, 2}, 
    tr_nominal::AbstractArray{D, 2},
    times::Vector{D}, 
    fieldmap::Vector{D},
    csm::Array{Complex{D}, 2};
    sim_method::BlochHighOrder=BlochHighOrder("111"), 
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOpv2w prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nSam, nCha)
    x0 = ones(Float64, nVox)
    if use_gpu
        out = out |> gpu
        tr_kspha = tr_kspha |> gpu
        tr_nominal = tr_nominal |> gpu
        xm = xm |> gpu
        x0 = x0 |> gpu
        x = grid.x |> gpu
        y = grid.y |> gpu
        z = grid.z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csm = csm |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(tr_kspha[1,p]), @view(tr_kspha[2,p]), @view(tr_kspha[3,p]), @view(tr_kspha[4,p]),
                                       @view(tr_kspha[5,p]), @view(tr_kspha[6,p]), @view(tr_kspha[7,p]), @view(tr_kspha[8,p]), @view(tr_kspha[9,p])
        hx, hy, hz = @view(tr_nominal[1,p]), @view(tr_nominal[2,p]), @view(tr_nominal[3,p])
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
    # out = out ./ sqrt(prod(grid.matrixSize))
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function ctprod_HighOrderOpv2w(
    ym::AbstractVector{T}, 
    grid::Grid{D},
    nVox::Int64, 
    nSam::Int64, 
    nCha::Int64,
    tr_kspha::AbstractArray{D, 2}, 
    tr_nominal::AbstractArray{D, 2},
    times::Vector{D}, 
    fieldmap::Vector{D},
    csm::Array{Complex{D}, 2};
    sim_method::BlochHighOrder=BlochHighOrder("111"), 
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    csmC = conj.(csm)
    ym = reshape(ym, nSam, nCha)
    # ym = Vector(ym[1:nSam])

    if verbose
        @info "HighOrderOpv2w ctprod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nVox, nCha)
    x0 = ones(Float64, nVox)

    if use_gpu
        out = out |> gpu
        tr_kspha = tr_kspha |> gpu
        tr_nominal = tr_nominal |> gpu
        ym = ym |> gpu
        x0 = x0 |> gpu
        x = grid.x |> gpu
        y = grid.y |> gpu
        z = grid.z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
        csmC = csmC |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(tr_kspha[1,p]), @view(tr_kspha[2,p]), @view(tr_kspha[3,p]), @view(tr_kspha[4,p]),
                                       @view(tr_kspha[5,p]), @view(tr_kspha[6,p]), @view(tr_kspha[7,p]), @view(tr_kspha[8,p]), @view(tr_kspha[9,p])
        hx, hy, hz = @view(tr_nominal[1,p]), @view(tr_nominal[2,p]), @view(tr_nominal[3,p])
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
    # out = out ./ sqrt(prod(grid.matrixSize))
    if use_gpu
        out = out |> cpu
    end
  return vec(sum(out, dims=2))
end

function Base.adjoint(op::HighOrderOpv2w{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
