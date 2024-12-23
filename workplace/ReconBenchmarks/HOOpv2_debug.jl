using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData
using BenchmarkTools
using CUDA
using MAT

CUDA.device!(1)

T = Float64
path         = "$(@__DIR__)/workplace/ReconBenchmarks"
# outpath   = "$(path)/out"; if ispath(outpath) == false mkpath(outpath) end

seq_file = "$(path)/data_r4_1p0/7T_1mm-200-r4_max51-fa90.seq"
mrd_file = "$(path)/data_r4_1p0/meas_MID00115_FID53003_pulseq_v0_r4_1p0_standard.mrd"
syn_file = "$(path)/data_r4_1p0/syn_meas_MID00115_FID53003_pulseq_v0_r4_1p0_standard.mat"
gre_file = "$(path)/data_r4_1p0/syn_meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mat"

# Coil-Sensitivity Map (CSM), ΔB0 map, mask
b0         = matread(gre_file)["b0"];
csm        = matread(gre_file)["csm"];    # (nY, nX, nCha)
mask       = matread(gre_file)["mask"];

# MRI signal data, FOV, matrix size, synchronized dynamic fields (kspha, rad, rad/m⁻¹, rad/m⁻²)
data       = matread(syn_file)["data"];
# data       = matread(syn_rep1_file)["data"]; 
FOV        = matread(syn_file)["FOV"];
matrixSize = matread(syn_file)["matrixSize"];
kStandard  = matread(syn_file)["ksphaStandard_syn"]./2π;
kStitched  = matread(syn_file)["ksphaStitched_syn"]./2π;

Δx, Δy, Δz = T.(FOV ./ matrixSize);
nX, nY, nZ = matrixSize;
nSample, nCha = size(data);
datatime = collect(1e-6 * (0:nSample-1));


# For HighOrderOp_v2
grid = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=true)
BHO = BlochHighOrder("111")
Nblocks = 5;
use_gpu = true;
verbose = true;

solver = "cgnr"; reg = "L2"; iter = 10; λ = 1e-3
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
weight = SampleDensity(kStitched'[2:3,:], (nX, nY));

########################################################################
# HighOrderOpv2
########################################################################
HOOpv2 = HighOrderOpv2(grid, T.(kStitched[:, 1:9]'), T.(datatime), BHO;
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time xv2 = recon_HOOp(HOOpv2, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(xv2))


HOOpv2_i2 = HighOrderOpv2_i2(grid, T.(-kStitched[:, 1:9]'), T.(datatime); sim_method=BHO,
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(-b0), use_gpu=use_gpu, verbose=verbose);  
@time xv2_i2 = recon_HOOp(HOOpv2_i2, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(xv2_i2))

kernel = [1/8 1/4 0 -1/4 -1/8]';
kStitched_dt = conv(kStitched,kernel)[3:end-2,:];
# plt_kspha(kStitched_dt, dt)
########################################################################
# HighOrderOpv2
########################################################################
grid        = grid
tr_kspha    = T.(-kStitched[:, 1:9]')
times       = T.(datatime)
fieldmap    = T.(-b0)
csm         = Complex{T}.(csm)
sim_method  = BlochHighOrder("111")
tr_nominal  = tr_kspha[2:4, :]
tr_kspha_dt = T.(-kStitched_dt[:, 1:9]')
Nblocks     = 5
use_gpu     = true
verbose     = false



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

@info "HighOrderOpv2_i2 Nblocks=$Nblocks, use_gpu=$use_gpu"
@info grid
@info sim_method


xm = vec(xv2)
ym = vec(data)


@time prod_HighOrderOpv2_i2(xm, bf, nVox, nSam, nCha, tr_kspha, times, fieldmap, csm;
    Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)

@time ctprod_HighOrderOpv2_i2(ym, bf, nVox, nSam, nCha, tr_kspha, times, fieldmap, csm;
    Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)


@benchmark prod_HighOrderOpv2_i2(xm, bf, nVox, nSam, nCha, tr_kspha, times, fieldmap, csm; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)
@benchmark ctprod_HighOrderOpv2_i2(ym, bf, nVox, nSam, nCha, tr_kspha, times, fieldmap, csm; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)
@benchmark prod_dt_HighOrderOpv2_i2(xm, bf, nVox, nSam, nCha, tr_kspha, tr_kspha_dt, times, fieldmap, csm; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu, verbose=verbose)


function prod_dt_HighOrderOpv2_i2(
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
        @info "HighOrderOpv2_i2 prod Nblocks=$Nblocks, use_gpu=$use_gpu"
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
    # out = out ./ sqrt(prod(grid.matrixSize))
    if use_gpu
        out = out |> cpu
    end
    return vec(out)
end