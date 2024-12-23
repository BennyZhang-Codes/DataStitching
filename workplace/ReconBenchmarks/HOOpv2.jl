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

# seq = read_seq(seq_file)[end-9:end-3]; 
# seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
# TE = seq.DEF["TE"];
# _, k_adc = get_kspace(seq);
# times = KomaMRIBase.get_adc_sampling_times(seq);
# times = times .- times[1] .+ TE;

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
verbose = false;

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
HOOpv2 = HighOrderOpv2(grid, T.(kStitched[:, 1:9]'), T.(datatime); sim_method=BHO,
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time xv2 = recon_HOOp(HOOpv2, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(xv2))

HOOpv2_i2 = HighOrderOpv2_i2(grid, T.(-kStitched[:, 1:9]'), T.(datatime); sim_method=BHO,
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(-b0), use_gpu=use_gpu, verbose=verbose);  
@time xv2_i2 = recon_HOOp(HOOpv2_i2, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(xv2_i2))


# Benckmark test
b1 = @benchmarkable recon_HOOp(   HOOpv2, Complex{T}.(data), Complex{T}.(weight), recParams) samples=10 evals=1 seconds=3600
t1 = run(b1)

b2 = @benchmarkable recon_HOOp(HOOpv2_i2, Complex{T}.(data), Complex{T}.(weight), recParams) samples=10 evals=1 seconds=3600
t2 = run(b2)

# Judgement
m1 = median(t1)
m2 = median(t2)
judgement = judge(m1, m2)







Op = HighOrderOp((nX, nY), tr_nominal, tr_Stitched, BlochHighOrder("111"); use_gpu=true, Δx=Δx, Δy=Δy, fieldmap=b0, Nblocks=Nblocks, grid=grid, verbose=verbose);

recParams[:encodingOps] = reshape([DiagOp(Op, nCha) ∘ S], 1,1);

@time rec = reconstruction(acqData, recParams).data[:,:];
plt_image(abs.(rec))
Nblocks=40;
Op_i2 = HighOrderOp_i2((nX, nY), tr_nominal, tr_Stitched, BlochHighOrder("011"); Δx=Δx, Δy=Δy, 
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=b0, grid=grid, use_gpu=true, verbose=verbose);

@time x = recon_HOOp(Op_i2, acqData, recParams)
plt_image(abs.(x))
plt_image(abs.(x-rec))



# Benckmark test
b1 = @benchmarkable recon_HOOp(Op_i2, acqData, recParams) samples=20 evals=1 seconds=3600
t1 = run(b1)

b2 = @benchmarkable reconstruction(acqData, recParams) samples=20 evals=1 seconds=36000
t2 = run(b2)

# Judgement
m1 = median(t1)
m2 = median(t2)
judgement = judge(m1, m2)

# # save benchmark test results
# matwrite("$(path)/data_r4_1p0/Benchmarktest_TITAN.mat", Dict("HOOp_i2" => t1, "HOOp" => t2, "judgement" => judgement))

# matwrite("$(path)/data_r4_1p0/Benchmarktest_3090.mat", Dict("HOOp_i2" => t1, "HOOp" => t2, "judgement" => judgement))

HOOpv2_i2 = HighOrderOpv2_i2(grid, T.(kStitched[:, 1:9]'), T.(datatime); sim_method=BHO,
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);





