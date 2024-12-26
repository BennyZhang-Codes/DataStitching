using KomaHighOrder
using MAT
using Interpolations
using PyPlot
using MRIReco
using RegularizedLeastSquares


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

HOOp = HighOrderOpv2_i2(gridding, T.(pha_spha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
# plt_image(abs.(x1))


# del0 = 0; # starting guess 
# maxNit_cgnr = length(resvec)-1;
# delJumpFact = 3;
# numCoarseSearch = 0;






# numCoarseSearch = 0;



# plt_kspha(kspha_dt, dt)
# plt_kspha(kspha, dt)
# plt_kspha_com(kspha_dt[1:end-1,:], diff(kspha, dims=1), dt)

# Get ready to start iterations
data            = data;
kspha           = ksphaStitched;
dt              = dt/dt
datatime        = datatime/1e-6
StartTime       = (delayStitched - 1e-6)/1e-6
csm             = csm;
b0              = b0 * 1e-6; 



τ = 0;
JumpFact = 3;
nIter = 1;
Δτ_prev = Δτ = Inf;
Δτ_min = 0.005 * dt; # us
@info "iterations..."
τ_perIter = Vector{Float64}();
push!(τ_perIter, τ);
verbose=false


using DSP
kernel = [1/8 1/4 0 -1/4 -1/8]';
kspha_dt = conv(kspha, kernel)[3:end-2,:]/dt;

while abs(Δτ) > Δτ_min
    # Interpolate to match datatimes
    kspha_τ    = InterpTrajTime(kspha   , dt, τ + StartTime)[1:nSample,1:9];
    kspha_dt_τ = InterpTrajTime(kspha_dt, dt, τ + StartTime)[1:nSample,1:9];
    
    HOOp    = HighOrderOpv2_i2(gridding, T.(kspha_τ'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    HOOp_dt = HighOrderOpv2_i2(gridding, T.(kspha_τ'), T.(datatime); sim_method=BHO, tr_kspha_dt=T.(kspha_dt_τ'), Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    
    # Update image
    x = recon_HOOp1(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams)
    # plt_image(abs.(x))
    
    # Update delay
    bhat = vec(data) - HOOp * vec(x); # Y - Aₚxₚ
    bhathat = HOOp_dt * vec(x);     # Bₚxₚ
    Δτ = JumpFact * real(bhathat \ bhat);
    τ += Δτ;
    @info "Iteration $nIter: Δτ = $(round(Δτ/dt, digits=5)) [us] | τ = $(round(τ, digits=5)) [us] | JumpFact = $JumpFact"
    if (nIter>1) && (sign(Δτ_prev) != sign(Δτ))
        println("Jumping back to previous solution")
        JumpFact = maximum([1, JumpFact/2]);
    end
    Δτ_prev = Δτ;
    nIter  += 1;
    push!(τ_perIter, τ);
end




function recon_HOOp1(HOOp::HighOrderOpv2_i2{Complex{T}}, Data::AbstractArray{Complex{T},2}, weight::AbstractVector{Complex{T}}, recParams::Dict) where T<:AbstractFloat
    recoParams = merge(defaultRecoParams(), recParams)
    reg = Regularization(recParams[:regularization], recParams[:λ]; shape=recParams[:reconSize])

    nSample, nCha = size(Data)
    Data = vec(Data)
    E = HOOp
    EᴴE = normalOperator(E)
    solver = createLinearSolver(recParams[:solver], E; AᴴA=EᴴE, reg=reg, recoParams...)
    x = solve(solver, Data; recoParams...)
    x = reshape(x, recParams[:reconSize])
    return x
end