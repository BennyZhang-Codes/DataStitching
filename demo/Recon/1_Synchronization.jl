using HighOrderMRI
using MAT, Interpolations
using CUDA
CUDA.device!(0)

T             = Float64;
path          = joinpath(@__DIR__, "demo/Recon")
data_file     = "$(path)/data_1p0_200_r4.mat"

@info "data file: $(data_file)"
data          = matread(data_file);

csm           = data["gre_csm"];            # coil sensitivity map
b0            = data["gre_b0"];             # ΔB0 map
mask          = data["gre_mask"];           # mask

kdata         = data["kdata"];              # spiral k-space data
datatime      = data["datatime"];           # time stamps of k-space data
matrixSize    = data["matrixSize"];         # matrix size
FOV           = data["FOV"];                # field of view

k0_ecc        = data["k0_ecc"];             # b0 compensation of scanner from ECC model
dt            = data["dfc_dt"];             # dwell time of field dynamics
ksphaStitched = data["dfc_ksphaStitched"];  # coefficients of the field dynamics with data stitching method
ksphaStandard = data["dfc_ksphaStandard"];  # coefficients of the field dynamics with standard method
startStitched = data["dfc_startStitched"];  # start time of the field dynamics with data stitching method
startStandard = data["dfc_startStandard"];  # start time of the field dynamics with standard method

ksphaNominal  = data["ksphaNominal"];       # nominal kspace trajectory, kx, ky
startNominal  = data["startNominal"];       # start time of the nominal trajectory

# settings for synchronization between the MRI data and the field dynamics
Δx, Δy, Δz    = T.(FOV ./ matrixSize);
nX, nY, nZ    = matrixSize;
gridding      = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)

intermode     = AkimaMonotonicInterpolation()
iter_max      = 20;
JumpFact      = 6;
Δτ_min        = T.(0.001);
λ             = T.(0);

recon_terms   = "111";
nBlock        = 20;
use_gpu       = true;
verbose       = false;

kdata         = kdata ./ exp.(1im.*k0_ecc)';

tau_Stitched  = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStitched[:,:]), T.(startStitched/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tau_Stitched: $(tau_Stitched) us"

tau_Standard  = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStandard[:,:]), T.(startStandard/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tau_Standard: $(tau_Standard) us"

tau_Nominal   = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaNominal[:,:]),  T.(startNominal/dt),  T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), recon_terms=recon_terms, nBlock=nBlock, use_gpu=use_gpu, verbose=verbose)
@info "tau_Nominal: $(tau_Nominal) us"
