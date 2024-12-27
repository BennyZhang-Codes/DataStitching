export HighOrderOp
# recon
using LinearOperators

"""
    A Julia implementation of the expanded signal encoding model.
- This implementation using GPU with CUDA.jl to accelerate the calculation.
- If the GPU memory is not enough, the calculation can be divided into blocks.
"""
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
    HighOrderOp(grid::Grid{T}, tr_kspha::AbstractArray{T, 2}, times::AbstractVector{T}; fieldmap::Matrix{T} = zeros(T,(grid.nX, grid.nY)), csm::Array{Complex{T}, 3} = ones(Complex{T},(grid.nX, grid.nY)..., 1), sim_method::BlochHighOrder = BlochHighOrder("111"), tr_nominal::AbstractArray{T, 2} = tr_kspha[2:4, :], tr_kspha_dt = nothing, Nblocks::Int64 = 50, use_gpu::Bool = true, verbose::Bool = false)

# Description
    generates a `HighOrderOp` which explicitely evaluates the MRI Fourier HighOrder encoding operator.

# Arguments:
* `grid::Grid{T}`                   - grid object.
* `tr_kspha::AbstractArray{T, 2}`   - trajectory of k-space phase encoding.
* `times::AbstractVector{T}`        - time points for trajectory.

# Keywords:
* `fieldmap::Matrix{T}`             - fieldmap for off-resonance correction.
* `csm::Array{Complex{T}, 3}`       - coil sensitivity map.
* `sim_method::BlochHighOrder`      - BlochHighOrder simulation method.
* `tr_nominal::AbstractArray{T, 2}` - nominal trajectory.
* `tr_kspha_dt`                     - trajectory of k-space phase encoding for time-derivative.
* `Nblocks::Int64`                  - split trajectory into `Nblocks` blocks to avoid memory overflow.
* `use_gpu::Bool`                   - use GPU for HighOrder encoding/decoding(default: `true`).
* `verbose::Bool`                   - print progress information(default: `false`).
"""
function HighOrderOp(
    grid        ::Grid{T},
    tr_kspha    ::AbstractArray{T, 2}, 
    times       ::AbstractVector{T};
    fieldmap    ::AbstractArray{T, 2}  = zeros(T,(grid.nX, grid.nY)), 
    csm         ::Array{Complex{T}, 3} = ones(Complex{T},(grid.nX, grid.nY)..., 1), 
    sim_method  ::BlochHighOrder       = BlochHighOrder("111"),
    tr_nominal  ::AbstractArray{T, 2}  = tr_kspha[2:4, :],
    tr_kspha_dt                        = nothing,
    Nblocks     ::Int64                = 50, 
    use_gpu     ::Bool                 = true, 
    verbose     ::Bool                 = false, 
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
    tr_kspha = prep_kspha(tr_kspha, tr_nominal, sim_method)

    if use_gpu
        tr_kspha    = tr_kspha    |> gpu
        tr_kspha_dt = tr_kspha_dt |> gpu
        grid        = grid        |> gpu
        times       = times       |> gpu
        fieldmap    = fieldmap    |> gpu
        csm         = csm         |> gpu
    end
    bf = basisfunc(grid.x, grid.y, grid.z, collect(1:9))
    
    if verbose
        @info "HighOrderOp Nblocks=$Nblocks, use_gpu=$use_gpu"
    end

    if isnothing(tr_kspha_dt)
        func_prod = (res,xm)->(res .= prod_HighOrderOp(xm, bf, nVox, nSam, nCha, tr_kspha, times, fieldmap, csm; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
    else
        @assert size(tr_kspha_dt) == size(tr_kspha) "tr_kspha_dt must have same size as tr_kspha"
        func_prod = (res,xm)->(res .= prod_dt_HighOrderOp(xm, bf, nVox, nSam, nCha, tr_kspha, tr_kspha_dt, times, fieldmap, csm; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
    end
    func_ctprod = (res,ym)->(res .= ctprod_HighOrderOp(ym, bf, nVox, nSam, nCha, tr_kspha, times, fieldmap, csm; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
    
    return HighOrderOp{Complex{T},Nothing,Function}(nRow, nCol, false, false
                , func_prod
                , nothing
                , func_ctprod
                , 0,0,0, false, false, false, Complex{T}[], Complex{T}[])
end

function prod_dt_HighOrderOp(
    xm::AbstractVector{T}, 
    bf::AbstractArray{D, 2},
    nVox::Int64,
    nSam::Int64, 
    nCha::Int64,
    tr_kspha::AbstractArray{D, 2}, 
    tr_kspha_dt::AbstractArray{D, 2},
    times::AbstractVector{D}, 
    fieldmap::AbstractVector{D},
    csm::AbstractArray{Complex{D}, 2};
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp prod_dt Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    if use_gpu
        xm = xm |> gpu
        out = CUDA.zeros(Complex{D}, nSam, nCha)
    else
        out = zeros(Complex{D}, nSam, nCha)
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        ϕ = @view(times[p]) .* fieldmap' .+ (bf * @view(tr_kspha[:,p]))'
        e = exp.(2*1im*pi*ϕ) .* (bf * @view(tr_kspha_dt[:,p]))' .* (2*1im*pi)
        out[p, :] =  e * (xm .* csm)
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
    out = out ./ sqrt(nVox)
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function prep_kspha(
    tr_kspha::AbstractArray{T, 2}, 
    tr_nominal::AbstractArray{T, 2}, 
    sim_method::BlochHighOrder) where T<:AbstractFloat
    if sim_method.ho0 == false
        tr_kspha[1, :] = tr_kspha[1, :] .* 0
    end
    if sim_method.ho1 == false
        tr_kspha[2:4, :] = tr_nominal[:, :]
    end
    if sim_method.ho2 == false
        tr_kspha[5:9, :] = tr_kspha[5:9, :] .* 0
    end
    return tr_kspha
end

function prod_HighOrderOp(
    xm::AbstractVector{T}, 
    bf::AbstractArray{D, 2},
    nVox::Int64,
    nSam::Int64, 
    nCha::Int64,
    tr_kspha::AbstractArray{D, 2}, 
    times::AbstractVector{D}, 
    fieldmap::AbstractVector{D},
    csm::AbstractArray{Complex{D}, 2};
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    if use_gpu
        xm = xm |> gpu
        out = CUDA.zeros(Complex{D}, nSam, nCha)
    else
        out = zeros(Complex{D}, nSam, nCha)
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        ϕ = @view(times[p]) .* fieldmap' .+ (bf * @view(tr_kspha[:,p]))'
        e = exp.(2*1im*pi*ϕ)
        out[p, :] =  e * (xm .* csm)
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end

    out = out ./ sqrt(nVox)
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function ctprod_HighOrderOp(
    ym::AbstractVector{T}, 
    bf::AbstractArray{D, 2},
    nVox::Int64, 
    nSam::Int64, 
    nCha::Int64,
    tr_kspha::AbstractArray{D, 2}, 
    times::AbstractVector{D}, 
    fieldmap::AbstractVector{D},
    csm::AbstractArray{Complex{D}, 2};
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nSam], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    csmC = conj.(csm)
    ym = reshape(ym, nSam, nCha)
    # ym = Vector(ym[1:nSam])

    if verbose
        @info "HighOrderOp ctprod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    
    if use_gpu
        ym = ym |> gpu
        out = CUDA.zeros(Complex{D}, nVox, nCha)
    else
        out = zeros(Complex{D}, nVox, nCha)
    end

    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        ϕ =  fieldmap .* @view(times[p])' .+ (bf * @view(tr_kspha[:,p]))
        e = exp.(2*1im*pi*ϕ)
        out +=  conj(e) * ym[p, :]
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end

    out = out ./ sqrt(nVox)
    out = out .* csmC
    if use_gpu
        out = out |> cpu
    end
  return vec(sum(out, dims=2))
end

function Base.adjoint(op::HighOrderOp{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
