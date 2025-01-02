using MAT, Interpolations


path         = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/"
out_path     = "$(path)/data"
if ispath(out_path) == false mkpath(out_path) end

seqs = ["7T_1mm-200-r4_max51-fa90.seq",
        "7T_1mm-200-r3_max51-fa90.seq",
        "7T_1mm-200-r2_max51-fa90.seq",
        "7T_0.71mm-280-r2_max51-fa90.seq",
        "7T_0.6mm-332-r3_max51-fa90.seq",
        "7T_0.5mm-400-r4_max51-fa90.seq"]

datas = ["1p0_200_r4.mat",
         "1p0_200_r3.mat",
         "1p0_200_r2.mat",
         "0p71_280_r2.mat",
         "0p6_332_r3.mat",
         "0p5_400_r4.mat"]

T = Float64;


idx = 1
seq = seqs[idx]
data = datas[idx]

seq_file = "$(path)/seq/$(seq)" 
data_file = "$(path)/data/$(data)"
@info "seq file: $(seq_file)"
@info "data file: $(dfc_file)"
seq  = read_seq(seq_file);
data = matread(data_file);

    
csm           = data["gre_csm"];
b0            = data["gre_b0"];
mask          = data["gre_mask"];

kdata         = data["kdata"];
matrixSize    = data["matrixSize"];
FOV           = data["FOV"];

dt            = data["dfc_dt"];
ksphaStitched = data["dfc_ksphaStitched"]./2π;  # rad to Hz
ksphaStandard = data["dfc_ksphaStandard"]./2π;
startStitched = data["dfc_startStitched"];
startStandard = data["dfc_startStandard"];

# plt_kspha_com(ksphaStitched, ksphaStandard, dt)
# nSample, nCha = size(kdata);
# pha_spha = InterpTrajTime(ksphaStitched,dt,startStitched)[1:nSample,:];
# datatime = collect(dt * (0:nSample-1));

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
Nblocks   = 10;
use_gpu   = true;
verbose   = false;

τStitched = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStitched[:,1:9]), T.(startStitched/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)

τStandard = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStandard), T.(startStandard/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)


import KomaHighOrder: get_grads
for s in seq
    if KomaHighOrder.is_ADC_on(s)
        t = collect(0:dt:s.DUR);
        gx, gy, gz = get_grads(s, t)
        kx = cumtrapz(ones(length(t)-1)'*dt, gx')
        ky = cumtrapz(ones(length(t)-1)'*dt, gy')
        kz = cumtrapz(ones(length(t)-1)'*dt, gz')
        k = [kx' ky' kz'] * γ
    end
end

τNominal = FindDelay_v2(gridding, Complex{T}.(kdata), T.(ksphaStandard), T.(startStandard/dt), T.(dt/dt);
            intermode=AkimaMonotonicInterpolation(), JumpFact=JumpFact, iter_max=iter_max, Δτ_min=Δτ_min, λ=λ,
            fieldmap=T.(b0*dt), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)
