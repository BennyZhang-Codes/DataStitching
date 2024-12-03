#=
only one calculation of the exponentials of total phase is needed, which is the most expensive part of the computation.
but it's tow slow with frequently tansfering data between CPU and GPU.
=#

export HighOrderOp_i1
# recon
using LinearOperators

mutable struct HighOrderOp_i1{T,F1,F2} <: AbstractLinearOperator{T}
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

LinearOperators.storage_type(op::HighOrderOp_i1) = typeof(op.Mv5)

"""
    HighOrderOp_i1(shape::NTuple{D,Int64}, tr_nominal::Trajectory, tr_measured::Trajectory; Nblocks::Int64=50, use_gpu::Bool=true) where D

# Description
    generates a `HighOrderOp_i1` which explicitely evaluates the MRI Fourier HighOrder encoding operator.

# Arguments:
* `shape::NTuple{D,Int64}`  - size of image to encode/reconstruct
* `tr_nominal::Trajectory`  - (must be 3 rows [x, y, z]') nominal Trajectory without normalization.
* `tr_measured::Trajectory`    - (must be 9 rows [h0, h1, h2, h3, h4, h5, h6, h7, h8]') measured Trajectory without normalization. 
``
                              this results in complex valued image even for real-valued input.
* `Nblocks`                 - split trajectory into `Nblocks` blocks to avoid memory overflow.
* `use_gpu`                 - use GPU for HighOrder encoding/decoding(default: `true`).
"""
function HighOrderOp_i1(
    shape::NTuple{D,Int64}, 
    tr_nominal::Trajectory, 
    tr_measured::Trajectory, 
    sim_method::BlochHighOrder; 
    fieldmap::Matrix{Float64}=zeros(Float64,shape), 
    Nblocks::Int64=50, 
    use_gpu::Bool=true, 
    verbose::Bool=false, 
    Δx::Float64=1e-3, 
    Δy::Float64=1e-3,
    grid::Int64=1,
    ) where D

    nodes_measured = Float64.(kspaceNodes(tr_measured))
    nodes_nominal = Float64.(kspaceNodes(tr_nominal))
    times = Float64.(readoutTimes(tr_measured))
    @assert size(nodes_measured,1) == 9 "nodes for measured must have 9 rows"
    @assert size(nodes_nominal,1) == 3 "nodes for nominal must have 3 rows"
    @assert size(fieldmap) == shape "fieldmap must have same size as shape"

    fieldmap = vec(fieldmap)
    nrow = size(nodes_measured,2)
    ncol = prod(shape)
    k = size(nodes_measured,2) # number of nodes
    Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k
    
    n = k÷Nblocks # number of nodes per block
    parts = [n for i=1:Nblocks] # number of nodes per block
    parts = [1+n*(i-1):n*i for i=1:Nblocks]
    if k%Nblocks!= 0
        push!(parts, n*Nblocks+1:k)
    end

    Nx, Ny = shape
    x, y = 1:Nx, 1:Ny
    if grid == 1      # x up->down, y left->right
        x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1
    elseif grid == 2
        x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
    elseif grid == 3  # x left->right, y down->up
        x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
    elseif grid == 4  # x left->right, y up->down
        x, y, z = vec(x*0.0 .+ y'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
    elseif grid == 5  # x right->left, y up->down
        x, y, z = vec(x*0.0 .+ y[end:-1:1]'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
        x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
    end
    # print(x)
    x, y = x * Δx, y * Δy 

    e = zeros(ComplexF64, nrow, ncol)
    calc_phase(@view(e[:,:]), nrow, ncol, x, y, z, nodes_measured, nodes_nominal, times, fieldmap;
        sim_method, Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)

    @info "HighOrderOp_i1 Nblocks=$Nblocks, use_gpu=$use_gpu, Δx=$Δx, Δy=$Δy"
    @info sim_method
    return HighOrderOp_i1{ComplexF64,Nothing,Function}(nrow, ncol, false, false
                , (res,xm)->(res .= prod_HighOrderOp_i1(@view(e[:,:]), xm, nrow, ncol;
                                        Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
                , nothing
                , (res,ym)->(res .= ctprod_HighOrderOp_i1(@view(e[:,:]), ym, nrow, ncol;
                                        Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose))
                , 0,0,0, false, false, false, ComplexF64[], ComplexF64[])
end


function calc_phase(
    e::AbstractArray{T}, 
    nrow::Int64,
    ncol::Int64, 
    x::Vector{Float64}, 
    y::Vector{Float64}, 
    z::Vector{Float64}, 
    nodes_measured::Matrix{Float64}, 
    nodes_nominal::Matrix{Float64},
    times::Vector{Float64}, 
    fieldmap::Vector{Float64};
    sim_method::BlochHighOrder=BlochHighOrder("111"), 
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:size(nodes_measured,2)], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where T<:Union{Real,Complex}

    @info "HighOrderOp_li calculaiton of the exponentials of total phase, prod Nblocks=$Nblocks, use_gpu=$use_gpu"

    x0 = ones(Float64, ncol)
    if use_gpu
        println("using GPU")
        nodes_measured = nodes_measured |> gpu
        nodes_nominal = nodes_nominal |> gpu
        x0 = x0 |> gpu
        x = x |> gpu
        y = y |> gpu
        z = z |> gpu
        times = times |> gpu
        fieldmap = fieldmap |> gpu
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
        @time exp_ϕ = exp.(-2*1im*pi*ϕ)
        if use_gpu
            exp_ϕ = exp_ϕ |> cpu
        end
        e[p, :] =  exp_ϕ
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
end


function prod_HighOrderOp_i1(
    e::AbstractArray{T}, 
    xm::AbstractVector{T}, 
    nrow::Int64,
    ncol::Int64;
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nrow], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where T<:Union{Real,Complex}
    xm = Vector(xm)
    if verbose
        @info "HighOrderOp_i1 prod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, nrow)
    if use_gpu
        out = out |> gpu
        xm = xm |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        exp_ϕ = e[p, :]
        if use_gpu
            exp_ϕ = exp_ϕ |> gpu
        end
        out[p] =  exp_ϕ * xm
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end

function ctprod_HighOrderOp_i1(
    e::AbstractArray{T}, 
    xm::AbstractVector{T}, 
    nrow::Int64,
    ncol::Int64;
    Nblocks::Int64=1, 
    parts::Vector{UnitRange{Int64}}=[1:nrow], 
    use_gpu::Bool=false, 
    verbose::Bool=false) where T<:Union{Real,Complex}

    xm = Vector(xm)
    if verbose
        @info "HighOrderOp_i1 ctprod Nblocks=$Nblocks, use_gpu=$use_gpu"
    end
    out = zeros(ComplexF64, ncol)

    if use_gpu
        out = out |> gpu
        xm = xm |> gpu
    end
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        exp_ϕ = transpose(e[p, :])
        if use_gpu
            exp_ϕ = exp_ϕ |> gpu
        end
        out +=  conj(exp_ϕ) * xm[p]
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
    if use_gpu
        out = out |> cpu
    end
  return vec(out)
end

function Base.adjoint(op::HighOrderOp_i1{T}) where T
  return LinearOperator{T}(op.ncol, op.nrow, op.symmetric, op.hermitian,
                        op.ctprod!, nothing, op.prod!)
end
