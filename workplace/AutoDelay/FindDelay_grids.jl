using KomaHighOrder
using MAT
using Interpolations
using PyPlot
using MRIReco
using RegularizedLeastSquares
import DSP: conv
using CUDA

CUDA.device!(7)

T=Float64;
path    = "$(@__DIR__)/workplace/AutoDelay/data"
outpath = "$(@__DIR__)/workplace/AutoDelay/out"

grefile = "$(path)/syn_meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mat"
MRIfile = "$(path)/syn_meas_MID00115_FID53003_pulseq_v0_r4_1p0_standard.mat"
DFCfile = "$(path)/7T_1p0_200_r4.mat"


csm  = matread(grefile)["csm"];
mask = matread(grefile)["mask"];
b0   = -matread(grefile)["b0"];

data       = matread(MRIfile)["data"];
matrixSize = matread(MRIfile)["matrixSize"];
FOV        = matread(MRIfile)["FOV"];

dt            = matread(DFCfile)["dt"];
ksphaStitched = matread(DFCfile)["ksphaStitched"]./2π;
ksphaStandard = matread(DFCfile)["ksphaStandard"]./2π;
delayStitched = matread(DFCfile)["delayStitched"];
delayStandard = matread(DFCfile)["delayStandard"];

plt_kspha_com(ksphaStitched, ksphaStandard, dt)

Δx, Δy, Δz = T.(FOV ./ matrixSize);
nX, nY, nZ = matrixSize;
nSample, nCha = size(data);
pha_spha = InterpTrajTime(ksphaStitched,dt,delayStitched)[1:nSample,:];
datatime = collect(dt * (0:nSample-1));


# For HighOrderOp_v2
gridding = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)
BHO = BlochHighOrder("111")
Nblocks = 10;
use_gpu = true;
verbose = false;

solver = "cgnr"; reg = "L2"; iter = 20; λ = 0.
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

weight = SampleDensity(pha_spha'[2:3,:], (nX, nY));

# HOOp = HighOrderOp(gridding, T.(pha_spha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
# plt_image(abs.(x1), vmaxp=99.9)

# HOOp = HighOrderOp(gridding, T.(pha_spha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp, Complex{T}.(data), recParams);
# plt_image(abs.(x1), vmaxp=99.9)

τ = FindDelay(gridding, Complex{T}.(data), T.(ksphaStitched), (delayStitched-dt)/dt, dt/dt;
    intermode=AkimaMonotonicInterpolation(), JumpFact=3, iter_max=20,
    fieldmap=T.(b0)*dt, csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)

τ = FindDelay_v2(gridding, Complex{T}.(data), T.(ksphaStitched), (delayStitched-dt)/dt, dt/dt;
    intermode=AkimaMonotonicInterpolation(), JumpFact=3, iter_max=20,
    fieldmap=T.(b0)*dt, csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)
    
######################################################################################
# Simulation with grids to visualize the delay effect
######################################################################################
# csm = csm_Birdcage(nX, nY, nCha);


im_grid = zeros(nX, nY);
for i in 10:10:(nX-1) # 绘制水平网格线
    im_grid[i, :] .= 1
end
for j in 10:10:(nY-1) # 绘制垂直网格线
    im_grid[:, j] .= 1
end
fig = plt_image(rotr90(abs.(im_grid)), vmaxp=99.9)
fig.savefig("$(outpath)/im_grid_200x200_7T_1p0_r4/fig_im_grid.png", dpi=600, transparent=false, bbox_inches="tight", pad_inches=0)
# plt_image(im_grid)


kspha = InterpTrajTime(ksphaStitched, dt, delayStitched, datatime);
HOOp = HighOrderOp(gridding, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
data = HOOp * vec(im_grid);
im1 = HOOp' * data;
data = reshape(data, nSample, nCha);
signalAmpl = sum(abs.(data), dims=1)/ nSample;


for snr in [Inf, 10, 5]
    for delay_applied in [-5.0, -2.5, -0.5, 0, 0.5, 2.5, 5.0];
        if snr == Inf
            data_in = data
        else
            data_in = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));
        end
        # recon with delay applied trajectories
        kspha = InterpTrajTime(ksphaStitched, dt, delayStitched - delay_applied*dt, datatime);
        weight = SampleDensity(kspha[:,2:3]', (nX, nY));
        HOOp = HighOrderOp(gridding, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
        @time x1 = recon_HOOp(HOOp, Complex{T}.(data_in), Complex{T}.(weight), recParams);
        fig = plt_image(abs.(x1), vmin=0, vmax=1)
        fig.savefig("$(outpath)/im_grid_200x200_7T_1p0_r4/fig_im_grid_snr$(snr)_del_$(round(delay_applied, digits=2)).png", dpi=600, transparent=false, bbox_inches="tight", pad_inches=0)
    end
end
