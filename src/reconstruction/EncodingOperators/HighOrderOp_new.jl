export HighOrderOp

using LinearOperators

"""
    A Julia implementation of the expanded signal encoding model.
- This implementation using GPU with CUDA.jl to accelerate the calculation.
- If the GPU memory is not enough, the calculation can be divided into blocks.
"""
mutable struct HighOrderOp{T,F1,F2} <: AbstractLinearOperator{T}
  nrow       :: Int
  ncol       :: Int
  symmetric  :: Bool
  hermitian  :: Bool
  prod!      :: Function
  tprod!     :: F1
  ctprod!    :: F2
  nprod      :: Int
  ntprod     :: Int
  nctprod    :: Int
  args5      :: Bool
  use_prod5! :: Bool
  allocated5 :: Bool
  Mv5        :: Vector{T}
  Mtu5       :: Vector{T}
end
LinearOperators.storage_type(op::HighOrderOp) = typeof(op.Mv5)


"""
    HighOrderOp(grid::Grid{T}, kspha::AbstractArray{T, 2}, times::AbstractVector{T}; fieldmap::Matrix{T} = zeros(T,(grid.nX, grid.nY)), csm::Array{Complex{T}, 3} = ones(Complex{T},(grid.nX, grid.nY)..., 1), sim_method::BlochHighOrder = BlochHighOrder("111"), k_nominal::AbstractArray{T, 2} = kspha[2:4, :], kspha_dt = nothing, nBlock::Int64 = 50, use_gpu::Bool = true, verbose::Bool = false)

# Description
    generates a `HighOrderOp` which explicitely evaluates the MRI Fourier HighOrder encoding operator.

# Arguments:
* `grid::Grid{T}`                   - grid object.
* `kspha::AbstractArray{T, 2}`      - [nSam, nTerm], Coefficients of field dynamics.
* `times::AbstractVector{T}`        - [nSam], time points for trajectory.

# Keywords:
* `fieldmap::Matrix{T}`             - [nX, nY], fieldmap for off-resonance correction.
* `csm::Array{Complex{T}, 3}`       - [nX, nY, nCha], coil sensitivity map.
* `recon_terms::BlochHighOrder`     - terms to be used in the Op.
* `k_nominal::AbstractArray{T, 2}`  - [nSam, 3], nominal kspace trajectory.
* `kspha_dt`                        - [nSam, nTerm], Coefficients of field dynamics for time-derivative.
* `nBlock::Int64`                   - split trajectory into `nBlock` blocks to avoid memory overflow.
* `use_gpu::Bool`                   - use GPU for HighOrder encoding/decoding(default: `true`).
* `verbose::Bool`                   - print progress information(default: `false`).
"""
function HighOrderOp(
    grid        :: Grid{T}                                                          ,
    kspha       :: AbstractArray{T, 2}                                              , 
    kcoco       :: AbstractArray{T, 2}                                              ,
    times       :: AbstractVector{T}                                                ;
    fieldmap    :: AbstractArray{T, 2}  = zeros(T,(grid.nX, grid.nY))               , 
    csm         :: Array{Complex{T}, 3} = ones(Complex{T},(grid.nX, grid.nY)..., 1) , 
    sim_method  :: BlochHighOrder       = BlochHighOrder("111")                     ,
    k_nominal   :: AbstractArray{T, 2}  = kspha[2:4, :]                             ,
    kspha_dt                            = nothing                                   ,
    nBlock      :: Int64                = 50                                        , 
    use_gpu     :: Bool                 = true                                      , 
    verbose     :: Bool                 = false                                     , 
    ) where {T<:AbstractFloat}
    nX, nY, nZ = grid.nX, grid.nY, grid.nZ
    nTerm, nSample = size(kspha)
    @assert nTerm in [9, 16] "kspha must have 9 or 16 terms (row) for up to 2nd or 3rd order terms"
    @assert size(k_nominal,1) == 3 "k_nominal must have 3 terms (row) for kx, ky, kz"
    @assert size(fieldmap) == (nX, nY) "FieldMap must have same size as $((nX, nY)) in grid"
    @assert size(csm)[1:2] == (nX, nY) "Coil-SensitivityMap must have same size as $((nX, nY)) in grid"
    
    nCha = size(csm, 3)
    nRow = nSample * nCha
    nCol = nVox = prod(grid.matrixSize)
    
    csm      = reshape(csm, nX*nY, nCha)   # [nX * nY, nCha]
    fieldmap = vec(fieldmap)

    k = nSam                               # number of samples
    nBlock = nBlock > k ? k : nBlock       # nBlock must be <= k
    
    n = k÷nBlock                           # number of sampless per block
    parts = [n for i=1:nBlock]             # number of samples per block
    parts = [1+n*(i-1):n*i for i=1:nBlock]
    if k%nBlock!= 0
        push!(parts, n*nBlock+1:k)
    end
    kspha = prep_kspha(kspha, k_nominal, sim_method)

    if use_gpu
        kspha       = kspha       |> gpu
        kspha_dt    = kspha_dt    |> gpu
        grid        = grid        |> gpu
        times       = times       |> gpu
        fieldmap    = fieldmap    |> gpu
        csm         = csm         |> gpu
    end
    bf = basisfunc(grid.x, grid.y, grid.z, collect(1:nTerm))
    
    if verbose
        @info "HighOrderOp nBlock=$nBlock, use_gpu=$use_gpu"
    end

    if isnothing(kspha_dt)
        func_prod = (res,xm)->(res .= prod_HighOrderOp(xm, bf, nVox, nSam, nCha, kspha, times, fieldmap, csm; nBlock=nBlock, parts=parts, use_gpu=use_gpu, verbose=verbose))
    else
        @assert size(kspha_dt) == size(kspha) "kspha_dt must have same size as kspha"
        func_prod = (res,xm)->(res .= prod_dt_HighOrderOp(xm, bf, nVox, nSam, nCha, kspha, kspha_dt, times, fieldmap, csm; nBlock=nBlock, parts=parts, use_gpu=use_gpu, verbose=verbose))
    end
    func_ctprod = (res,ym)->(res .= ctprod_HighOrderOp(ym, bf, nVox, nSam, nCha, kspha, times, fieldmap, csm; nBlock=nBlock, parts=parts, use_gpu=use_gpu, verbose=verbose))
    
    return HighOrderOp{Complex{T},Nothing,Function}(
                        nRow, nCol, 
                        false, false,
                        func_prod, nothing, func_ctprod,
                        0, 0, 0, 
                        false, false, false, 
                        Complex{T}[], Complex{T}[])
end


"""
    prod_dt_HighOrderOp
    for calculation of Bx (2022, https://doi.org/10.1002/mrm.29460)
"""
function prod_dt_HighOrderOp(
    xm        :: AbstractVector{T}                   , 
    bf        :: AbstractArray{D, 2}                 ,
    nVox      :: Int64                               ,
    nSam      :: Int64                               , 
    nCha      :: Int64                               ,
    kspha     :: AbstractArray{D, 2}                 , 
    kspha_dt  :: AbstractArray{D, 2}                 ,
    times     :: AbstractVector{D}                   , 
    fieldmap  :: AbstractVector{D}                   ,
    csm       :: AbstractArray{Complex{D}, 2}        ;
    nBlock    :: Int64                    = 1        , 
    parts     :: Vector{UnitRange{Int64}} = [1:nSam] , 
    use_gpu   :: Bool                     = false    , 
    verbose   :: Bool                     = false    ,
    ) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp prod_dt nBlock=$nBlock, use_gpu=$use_gpu"
    end
    if use_gpu
        xm = xm |> gpu
        out = CUDA.zeros(Complex{D}, nSam, nCha)
    else
        out = zeros(Complex{D}, nSam, nCha)
    end
    progress_bar = Progress(nBlock)
    for (block, p) = enumerate(parts)
        ϕ = @view(times[p]) .* fieldmap' .+ (bf * @view(kspha[:,p]))'
        e = exp.(2*1im*pi*ϕ) .* (bf * @view(kspha_dt[:,p]))' .* (2*1im*pi)
        out[p, :] =  e * (xm .* csm)
        if verbose
            next!(progress_bar, showvalues=[(:nBlock, block)])
        end
    end
    out = out ./ sqrt(nVox)
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function prep_kspha(
    kspha::AbstractArray{T, 2}, 
    k_nominal::AbstractArray{T, 2}, 
    sim_method::BlochHighOrder) where T<:AbstractFloat
    if sim_method.ho0 == false
        kspha[1, :] = kspha[1, :] .* 0
    end
    if sim_method.ho1 == false
        kspha[2:4, :] = k_nominal[:, :]
    end
    if sim_method.ho2 == false
        kspha[5:9, :] = kspha[5:9, :] .* 0
    end
    return kspha
end


"""
    Forward operator for HighOrderOp
"""
function prod_HighOrderOp(
    xm        :: AbstractVector{T}                   , 
    bf        :: AbstractArray{D, 2}                 ,
    nVox      :: Int64                               ,
    nSam      :: Int64                               , 
    nCha      :: Int64                               ,
    kspha     :: AbstractArray{D, 2}                 , 
    times     :: AbstractVector{D}                   , 
    fieldmap  :: AbstractVector{D}                   ,
    csm       :: AbstractArray{Complex{D}, 2}        ;
    nBlock    :: Int64                    = 1        , 
    parts     :: Vector{UnitRange{Int64}} = [1:nSam] , 
    use_gpu   :: Bool                     = false    , 
    verbose   :: Bool                     = false    ,
    ) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp prod nBlock=$nBlock, use_gpu=$use_gpu"
    end
    if use_gpu
        xm = xm |> gpu
        out = CUDA.zeros(Complex{D}, nSam, nCha)
    else
        out = zeros(Complex{D}, nSam, nCha)
    end
    progress_bar = Progress(nBlock)
    for (block, p) = enumerate(parts)
        ϕ = @view(times[p]) .* fieldmap' .+ (bf * @view(kspha[:,p]))'
        e = exp.(2*1im*pi*ϕ)
        out[p, :] =  e * (xm .* csm)
        if verbose
            next!(progress_bar, showvalues=[(:nBlock, block)])
        end
    end

    out = out ./ sqrt(nVox)
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end


"""
    Adjoint of prod_HighOrderOp
"""
function ctprod_HighOrderOp(
    ym        :: AbstractVector{T}                   , 
    bf        :: AbstractArray{D, 2}                 ,
    nVox      :: Int64                               , 
    nSam      :: Int64                               , 
    nCha      :: Int64                               ,
    kspha     :: AbstractArray{D, 2}                 , 
    times     :: AbstractVector{D}                   , 
    fieldmap  :: AbstractVector{D}                   ,
    csm       :: AbstractArray{Complex{D}, 2}        ;
    nBlock    :: Int64                    = 1        , 
    parts     :: Vector{UnitRange{Int64}} = [1:nSam] , 
    use_gpu   :: Bool                     = false    , 
    verbose   :: Bool                     = false    ,
    ) where {D<:AbstractFloat, T<:Union{Real,Complex}}
    csmC = conj.(csm)
    ym = reshape(ym, nSam, nCha)
    # ym = Vector(ym[1:nSam])

    if verbose
        @info "HighOrderOp ctprod nBlock=$nBlock, use_gpu=$use_gpu"
    end
    
    if use_gpu
        ym = ym |> gpu
        out = CUDA.zeros(Complex{D}, nVox, nCha)
    else
        out = zeros(Complex{D}, nVox, nCha)
    end

    progress_bar = Progress(nBlock)
    for (block, p) = enumerate(parts)
        ϕ =  fieldmap .* @view(times[p])' .+ (bf * @view(kspha[:,p]))
        e = exp.(2*1im*pi*ϕ)
        out +=  conj(e) * ym[p, :]
        if verbose
            next!(progress_bar, showvalues=[(:nBlock, block)])
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
  return LinearOperator{T}(
                            op.ncol, 
                            op.nrow, 
                            op.symmetric, 
                            op.hermitian,
                            op.ctprod!, 
                            nothing, 
                            op.prod!)
end
