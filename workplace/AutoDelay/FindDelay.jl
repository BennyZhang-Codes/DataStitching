using KomaHighOrder
using MAT
using Interpolations
using PyPlot
using MRIReco
using RegularizedLeastSquares
import DSP: conv

T=Float64;
path = "$(@__DIR__)/workplace/AutoDelay/data"

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
verbose = true;

solver = "cgnr"; reg = "L2"; iter = 10; λ = 0.
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

weight = SampleDensity(pha_spha'[2:3,:], (nX, nY));

# HOOp = HighOrderOp(gridding, T.(pha_spha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
# plt_image(abs.(x1))

τ = FindDelay1(gridding, Complex{T}.(data), T.(ksphaStitched), (delayStitched-dt)/dt, dt/dt;
    intermode=AkimaMonotonicInterpolation(), JumpFact=1,
    fieldmap=T.(b0)*dt, csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)

    
################################################################################
# debug 
################################################################################

gridding    = gridding                      
data        = Complex{T}.(data) 
kspha       = T.(ksphaStitched)           
StartTime   = delayStitched-dt                         
dt          = dt               
JumpFact    = 3      
intermode   = AkimaMonotonicInterpolation()
fieldmap    = T.(b0)
csm         = Complex{T}.(csm)
sim_method  = BlochHighOrder("111")
Nblocks     = 10     
use_gpu     = true   
solver      = "cgnr" 
reg         = "L2"   
iter_max    = 5     
λ           = 0.     
verbose     = true  



# 1. Initialize some variables
τ         = 0;
nIter     = 1;
Δτ_prev   = Δτ = Inf;
Δτ_min    = 0.005 * 1e-6;        # [us]
τ_perIter = Vector{Float64}();
push!(τ_perIter, τ);
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (gridding.nX, gridding.nY)
recParams[:regularization] = reg
recParams[:λ]              = λ
recParams[:iterations]     = iter_max
recParams[:solver]         = solver
datatime = T.(collect(dt * (0:nSample-1)));

# 2. Compute kspha_dt, which is "dk/dt"
kernel = [1/8 1/4 0 -1/4 -1/8]';
kspha_dt = conv(kspha, kernel)[3:end-2,:]/dt;

while abs(Δτ) > Δτ_min
    # Interpolate to match datatimes
    kspha_τ    = T.(InterpTrajTime(kspha   , dt, τ + StartTime, intermode=intermode)[1:nSample,1:9]');
    kspha_dt_τ = T.(InterpTrajTime(kspha_dt, dt, τ + StartTime, intermode=intermode)[1:nSample,1:9]');
    
    HOOp    = HighOrderOp(gridding, kspha_τ, datatime; sim_method=sim_method, 
                Nblocks=Nblocks, csm=csm, fieldmap=fieldmap, use_gpu=use_gpu, verbose=verbose);
    HOOp_dt = HighOrderOp(gridding, kspha_τ, datatime; sim_method=sim_method, tr_kspha_dt=kspha_dt_τ, 
                Nblocks=Nblocks, csm=csm, fieldmap=fieldmap, use_gpu=use_gpu, verbose=verbose);
    
    # Update image
    x = recon_HOOp(HOOp, data, recParams)
    if verbose
        plt_image(abs.(x), title="Iteration $nIter", vmaxp=99.9)
    end
    # Update delay
    y1 = vec(data) - HOOp * vec(x); # Y - Aₚxₚ
    y2 = HOOp_dt * vec(x);       # Bₚxₚ
    Δτ = JumpFact * real(y2 \ y1);
    τ += Δτ;
    @info "Iteration $nIter: Δτ = $(round(Δτ/dt, digits=5)) [us] | τ = $(round(τ/dt, digits=5)) [us] | JumpFact = $JumpFact"
    if (nIter>1) && (sign(Δτ_prev) != sign(Δτ))
        println("Jumping back to previous solution")
        JumpFact = maximum([1, JumpFact/2]);
    end
    Δτ_prev = Δτ;
    nIter  += 1;
    push!(τ_perIter, τ);
end

