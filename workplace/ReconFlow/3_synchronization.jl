using MAT, Interpolations
using CUDA
CUDA.device!(1)

T = Float64;

path         = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/"
out_path     = "$(path)/data"
if ispath(out_path) == false mkpath(out_path) end



datas = ["1p0_200_r4.mat",
         "1p0_200_r3.mat",
         "1p0_200_r2.mat",
         "0p71_280_r2.mat",
         "0p6_332_r3.mat",
         "0p5_400_r4.mat"]

for idx in [3, 2, 1]
idx = 6

data = datas[idx]
data_file = "$(path)/data/$(data)"
@info "data file: $(data_file)"
data = matread(data_file);

csm           = data["gre_csm"];
b0            = data["gre_b0"];
mask          = data["gre_mask"];

kdata         = data["kdata"];
datatime      = data["datatime"];
matrixSize    = data["matrixSize"];
FOV           = data["FOV"];

dt            = data["dfc_dt"];
ksphaStitched = data["dfc_ksphaStitched"];  
ksphaStandard = data["dfc_ksphaStandard"];
startStitched = data["dfc_startStitched"];
startStandard = data["dfc_startStandard"];

ksphaNominal  = data["ksphaNominal"];
startNominal  = data["startNominal"];

# FindDelay, using a model-based approach
Δx, Δy, Δz = T.(FOV ./ matrixSize);
nX, nY, nZ = matrixSize;

gridding  = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)

intermode = AkimaMonotonicInterpolation()
iter_max  = 20;
JumpFact  = 6;
Δτ_min    = T.(0.001);
λ         = T.(0);

BHO       = BlochHighOrder("111")
Nblocks   = 150;
use_gpu   = true;
verbose   = false;

tauStitched = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStitched[:,1:9]), T.(startStitched/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks, use_gpu=use_gpu, verbose=verbose)

tauStandard = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStandard[:,1:9]), T.(startStandard/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks, use_gpu=use_gpu, verbose=verbose)

tauNominal  = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaNominal[:,1:9]),  T.(startNominal/dt),  T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks, use_gpu=use_gpu, verbose=verbose)

# plt_kspha_com(ksphaNominal, InterpTrajTime(ksphaNominal, dt, τNominal*dt), dt)
# plt_kspha(ksphaNominal, dt)
# plt_kspha(ksphaStitched, dt)
# plt_kspha_com(ksphaNominal, ksphaStitched[1:29940,:], dt)


data["tauNominal"]      = tauNominal;
data["dfc_tauStitched"] = tauStitched;
data["dfc_tauStandard"] = tauStandard;


@info "Synchronization done" τ
MAT.matwrite(data_file, data)

end