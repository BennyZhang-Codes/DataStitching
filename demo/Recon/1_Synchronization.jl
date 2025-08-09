using HighOrderMRI
using MAT, Interpolations
using CUDA
CUDA.device!(0)

T             = Float64;
path          = joinpath(@__DIR__, "demo/Recon")
data_mat      = "7T_2D_Spiral_1p0_200_r4.mat"  
"""
7T_2D_Spiral_1p0_200_r4.mat
7T_2D_Spiral_0p5_400_r4.mat
7T_2D_EPI_1p0_200_r4.mat
7T_2D_EPI_0p5_400_r5.mat
""" 
data_file     = joinpath(path, data_mat)

@info "data file: $(data_file)"
data          = matread(data_file);

csm           = data["gre_csm"];            # coil sensitivity map
b0            = data["gre_b0"];             # ΔB0 map
mask          = data["gre_mask"];           # mask

kdata         = data["kdata"];              # spiral k-space data
datatime      = data["datatime"];           # time stamps of k-space data
matrixSize    = data["matrixSize"];         # matrix size
FOV           = data["FOV"];                # field of view

k0_ecc        = data["k0_adc"];             # b0 compensation of scanner from ECC model

dt_adc        = data["dt_adc"];

# dwell time of the trajectories
dt_Stitched   = data["dt_Stitched"];
dt_Standard   = data["dt_Standard"];
dt_GIRF      = data["dt_GIRF"];
dt_Nominal    = data["dt_Nominal"];

# phase coefficients
ksphaStitched = data["ksphaStitched"];     
ksphaStandard = data["ksphaStandard"];     
ksphaGIRF     = data["ksphaGIRF"];         
ksphaNominal  = data["ksphaNominal"];      

# start time of the field dynamics
startStitched = data["startStitched"]; 
startStandard = data["startStandard"]; 
startGIRF     = data["startGIRF"];     
startNominal  = data["startNominal"];  

datatime      = vec(datatime);

# settings for synchronization between the MRI data and the field dynamics
Δx, Δy, Δz    = T.(FOV ./ matrixSize);
nX, nY, nZ    = matrixSize;
gridding      = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)

intermode     = AkimaMonotonicInterpolation()
iter_max      = 20;
JumpFact      = 6;              # 6 for Spiral, 100 for EPI
Δτ_min        = T.(0.001);
λ             = T.(0);

nBlock        = 20;
use_gpu       = true;
verbose       = false;

dt = T.(dt_adc);
kdata = kdata ./ exp.(2π*1im.*k0_ecc)';

recon_terms = "1111"
ksphaStitched = InterpTrajTime(ksphaStitched, dt_Stitched, 0, collect(0:dt_adc:size(ksphaStitched, 1)*dt_Stitched));
tauStitched = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStitched[:,1:16]), T.(datatime/dt), T.(startStitched/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.((b0)*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tauStitched: $(tauStitched*dt*1e6) us"

recon_terms = "1111"
ksphaStandard = InterpTrajTime(ksphaStandard, dt_Standard, 0, collect(0:dt_adc:size(ksphaStandard, 1)*dt_Standard));
tauStandard = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStandard[:,1:16]), T.(datatime/dt), T.(startStandard/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.((b0)*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tauStandard: $(tauStandard*dt*1e6) us"

recon_terms = "1111"
ksphaGIRF  = InterpTrajTime(ksphaGIRF, dt_GIRF, 0, collect(0:dt_adc:size(ksphaGIRF, 1)*dt_GIRF));
tauGIRF    = FindDelay_v2(gridding, Complex{T}.(kdata), T.(     ksphaGIRF[:,1:16]), T.(datatime/dt), T.(   startGIRF/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.((b0)*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tauGIRF: $(tauGIRF*dt*1e6) us"

recon_terms = "010"  
ksphaNominal = InterpTrajTime(ksphaNominal, dt_Nominal, 0, collect(0:dt_adc:size(ksphaNominal, 1)*dt_Nominal));
tauNominal  = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaNominal[:,1:9]),  T.(datatime/dt),T.(startNominal/dt),  T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.((b0)*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tauNominal: $(tauNominal*dt*1e6) us"
