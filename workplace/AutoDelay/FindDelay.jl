using MAT
using Interpolations
using PyPlot
using MRIReco
using RegularizedLeastSquares


T=Float32;
path = "$(@__DIR__)/workplace/AutoDelay/data"

grefile = "$(path)/syn_meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mat"
MRIfile = "$(path)/syn_meas_MID00115_FID53003_pulseq_v0_r4_1p0_standard.mat"
DFCfile = "$(path)/7T_1p0_200_r4.mat"


csm  = matread(grefile)["csm"];
mask = matread(grefile)["mask"];
b0   = matread(grefile)["b0"];

data       = matread(MRIfile)["data"];
matrixSize = matread(MRIfile)["matrixSize"];
FOV        = matread(MRIfile)["FOV"];

dt            = matread(DFCfile)["dt"];
ksphaStitched = matread(DFCfile)["ksphaStitched"]./2π;
ksphaStandard = matread(DFCfile)["ksphaStandard"]./2π;
delayStitched = matread(DFCfile)["delayStitched"];
delayStandard = matread(DFCfile)["delayStandard"];

Δx, Δy, Δz = T.(FOV ./ matrixSize);
nX, nY, nZ = matrixSize;
nSample, nCha = size(data);
pha_spha = InterpTrajTime(ksphaStitched,dt,delayStitched-dt*0.5)[1:nSample,:];
datatime = collect(dt * (0:nSample-1));



# For HighOrderOp_v2
grid = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=true)
BHO = BlochHighOrder("010")
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

HOOp = HighOrderOpv2(grid, T.(pha_spha[:, 1:9]'), T.(datatime), BHO;
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);

weight = SampleDensity(pha_spha'[2:3,:], (nX, nY));
# x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams)
# plt_image(abs.(x1))


del0 = 0; # starting guess 
maxNit_cgnr = length(resvec)-1;
delJumpFact = 3;
numCoarseSearch = 0;




delSk0          = 0;
data_in         = kdata;
tdwell          = dt;
kspha           = InterpTrajTime(ksphaStitched,dt,delayStitched-dt);
datatime

b0              = b0; 
csm             = csm;
maxIter_cgnr    = 20;   
delJumpFact     = 3;
numCoarseSearch = 0;



kernel = [1/8 1/4 0 -1/4 -1/8]';
kspha_dt = conv(kspha,kernel)[3:end-2,:];
plt_kspha(kspha_dt, dt)
plt_kspha(kspha, dt)
plt_kspha_com(kspha_dt[1:end-1,:], diff(kspha, dims=1), dt)

# Get ready to start iterations
τ = delSk0;
JumpFact = 3;
nIter = 1;
Δτ = Inf;
Δτmin = 0.005; # us
@info "iterations..."
τ_perIter = zeros(1000,1);
τ_perIter[1] = τ;

while abs(Δτ) > Δτmin
    # Interpolate to match datatimes
    kspha_τ    = InterpTrajTime(kspha   , dt, τ)[1:nSample,1:9];
    kspha_dt_τ = InterpTrajTime(kspha_dt, dt, τ)[1:nSample,1:9];

    # Update image
    HOOp = HighOrderOpv2(grid, T.(kspha_τ'), T.(datatime), BHO; Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams)
    plt_image(abs.(x))

    HOOp = HighOrderOpv2(grid, T.(kspha_τ'), T.(datatime), BHO; tr_kspha_dt=T.(kspha_dt_τ'), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    
    # Update delay
    bhat = data_in - HOOp * x; # Y - Aₚxₚ
    bhathat                    # Bₚxₚ
    Δτ = JumpFact * real(bhathat \ bhat);
    τ += Δτ;
    @info "Iteration $nIter: Δτ = $Δτ, τ = $τ, JumpFact = $JumpFact"
    if (nIter>1) && (sign(Δ𝜏_prev) != sign(Δτ))
        JumpFact = maximum([1, JumpFact/2]);
    end
    Δτ_prev = Δτ;
    nIter += 1;
    τ_perIter[nIter] = τ;
end