using HighOrderMRI
using MAT
using PyPlot
using MRIReco
using BenchmarkTools
using CUDA

CUDA.device!(0)

T=Float64;
path = "$(@__DIR__)/workplace/AutoDelay/comparison_matmri"

data_file = "$(path)/data_images.mat"
traj_file = "$(path)/data_traj_spiral_R4.mat"
data_mat = matread(data_file)
traj_mat = matread(traj_file)

nSli = 1;
x    = data_mat["X"][:,:,nSli];
y    = data_mat["Y"][:,:,nSli];
z    = data_mat["Z"][:,:,nSli];
csm  = data_mat["Crcvr"][:,:,nSli,:];
b0   = data_mat["b0map"][:,:,nSli] ./ (2π);
im0  = data_mat["im0"][:,:,nSli];

datatime  = vec(traj_mat["datatime"]);
dt        = traj_mat["tdwelltraj"];
kspha_raw = traj_mat["phs_spha"][:, 1:9] ./ (2π);
StartTime = 50;

nX, nY, nCha = size(csm); nZ = 1
nSample = length(datatime)


##
g = Grid(nX=nX, nY=nY, nZ=nZ, Δx=1.5, Δy=1.5, Δz=1, x=vec(x), y=vec(y), z=vec(z))
BHO = BlochHighOrder("111")
Nblocks = 10;
use_gpu = true;
verbose = true;

kspha = InterpTrajTime(kspha_raw, dt, StartTime, datatime);
weight = SampleDensity(kspha[:,2:3]', (nX, nY));
HOOp_v3 = HighOrderOp_v3(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
data = HOOp_v3 * vec(im0);
im1 = HOOp_v3' * data;
data = reshape(data, nSample, nCha);


solver = "cgnr"; reg = "L2"; iter = 20; λ = 0.
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

# recon_terms = "010";
# HOOp_v3 = HighOrderOp_v3(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("010"), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp_v3, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x1)); title="HighOrderOp_v3", vmaxp=99.9)

# HOOp_new= HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); recon_terms=recon_terms, nBlock=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x2 = recon_HOOp(HOOp_new, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x2)); title="HighOrderOp", vmaxp=99.9)

# print(abs.(x1) ≈ abs.(x2))

# recon_terms = "110";
# HOOp_v3 = HighOrderOp_v3(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder(recon_terms), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp_v3, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x1)); title="HighOrderOp_v3", vmaxp=99.9)

# HOOp_new= HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); recon_terms=recon_terms, nBlock=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x2 = recon_HOOp(HOOp_new, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x2)); title="HighOrderOp", vmaxp=99.9)

# print(abs.(x1) ≈ abs.(x2))

# recon_terms = "111";
# HOOp_v3 = HighOrderOp_v3(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder(recon_terms), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp_v3, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x1)); title="HighOrderOp_v3", vmaxp=99.9)

# HOOp_new= HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); recon_terms=recon_terms, nBlock=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x2 = recon_HOOp(HOOp_new, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x2)); title="HighOrderOp", vmaxp=99.9)

# print(abs.(x1) ≈ abs.(x2))

Nblocks=10
verbose = true;
recon_terms = "111";
HOOp_v3 = HighOrderOp_v3(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder(recon_terms), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time recon_HOOp(HOOp_v3, Complex{T}.(data), Complex{T}.(weight), recParams);


HOOp= HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); recon_terms=recon_terms, nBlock=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams) ;


verbose = false;
recon_terms = "111";
HOOp_v3 = HighOrderOp_v3(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder(recon_terms), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
b1 = @benchmarkable recon_HOOp(HOOp_v3, Complex{T}.(data), Complex{T}.(weight), recParams) samples=20 evals=1 seconds=3600
t1 = run(b1)

HOOp_new= HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); recon_terms=recon_terms, nBlock=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
b2 = @benchmarkable recon_HOOp(HOOp_new, Complex{T}.(data), Complex{T}.(weight), recParams) samples=20 evals=1 seconds=3600
t2 = run(b2)

# Judgement
m1 = median(t1)
m2 = median(t2)
judgement = judge(m1, m2)


HOOp_new= HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); recon_terms=recon_terms, nBlock=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
b3 = @benchmarkable recon_HOOp(HOOp_new, Complex{T}.(data), Complex{T}.(weight), recParams) samples=20 evals=1 seconds=3600
t3 = run(b3)


k = rand(900000, 9)
x = vec(im0)
HOOp_v3 = HighOrderOp_v3(g, T.(k'), T.(collect(1:900000)); sim_method=BlochHighOrder(recon_terms), 
    csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=false, Nblocks=50);
@time HOOp_v3 * x;

HOOp_new = HighOrderOp(g, T.(k'), T.(collect(1:900000)); recon_terms=recon_terms,
    csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=false, nBlock=50);
@time HOOp_new * x;



nSam = 
nBlock = 50

nBlock = nBlock > nSam  ? nSam  : nBlock  # nBlock must be <= k
n      = nSam ÷ nBlock                    # number of sampless per block
parts  = [n for i=1:nBlock]               # number of samples per block
parts  = [1+n*(i-1):n*i for i=1:nBlock]
if nSam%nBlock!= 0
    push!(parts, n*nBlock+1:nSam)
end

nSam = 26
nBlock = 26

nBlock       = min(nBlock, nSam)                   # nBlock must be <= nSam
@assert nSam % nBlock == 0  "nSam should be evenly divided by nBlock."
nSamPerBlock = nSam ÷ nBlock
parts  = [1+nSamPerBlock*(i-1):nSamPerBlock*i for i=1:nBlock]