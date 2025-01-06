using MAT
# using CUDA
# CUDA.device!(7)

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

# for idx in [3, 2, 1]
idx = 1

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
ksphaNominal  = data["ksphaNominal"];

startStitched = data["dfc_startStitched"];
startStandard = data["dfc_startStandard"];
startNominal  = data["startNominal"];

tauStitched   = data["dfc_tauStitched"];
tauStandard   = data["dfc_tauStandard"];
tauNominal    = data["tauNominal"];
    
# plt_kspha_com(InterpTrajTime(ksphaStitched, dt, startStitched, datatime), InterpTrajTime(ksphaStitched, dt, startStitched+tauStitched*dt, datatime), dt)
# plt_kspha_com(InterpTrajTime(ksphaStandard, dt, startStandard, datatime), InterpTrajTime(ksphaStandard, dt, startStandard+tauStandard*dt, datatime), dt)
# plt_kspha_com(InterpTrajTime(ksphaNominal , dt, startNominal , datatime), InterpTrajTime(ksphaNominal , dt, startNominal +tauNominal *dt, datatime), dt)
ksphaStitched = InterpTrajTime(ksphaStitched, dt, startStitched+tauStitched*dt, datatime);
ksphaStandard = InterpTrajTime(ksphaStandard, dt, startStandard+tauStandard*dt, datatime);
ksphaNominal  = InterpTrajTime(ksphaNominal , dt, startNominal +tauNominal *dt, datatime);
# plt_kspha(ksphaStitched, dt)
# plt_kspha(ksphaStandard, dt)
# plt_kspha(ksphaNominal , dt)

# prepare some parameters for reconstruction
# 1. gridding 
Δx, Δy, Δz = T.(FOV ./ matrixSize);
nX, nY, nZ = matrixSize;
gridding  = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)

# 2. sampling density 
weightStitched = SampleDensity(ksphaStitched'[2:3,:], (nX, nY));
weightStandard = SampleDensity(ksphaStandard'[2:3,:], (nX, nY));
weightNominal  = SampleDensity( ksphaNominal'[2:3,:], (nX, nY));

use_gpu = true;
verbose = false;
Nblocks = 150;

solver = "cgnr"; reg = "L2"; iter = 20; λ = 1e-7
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

BHO = BlochHighOrder("111")
HOOp = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weightStitched), recParams);
plt_image(abs.(x), vmaxp=99.9)

# 1. stitched "111"
HOOp1 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("111"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 2. stitched "110"
HOOp2 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("110"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 3. stitched "011"
HOOp3 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("011"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 4. stitched "101"
HOOp4 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("101"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 5. stitched "100"
HOOp5 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("100"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 6. stitched "010"
HOOp6 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("010"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 7. stitched "001"
HOOp7 = HighOrderOp(gridding, T.(ksphaStitched[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("001"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);

HOOps          = [HOOp1, HOOp2, HOOp3, HOOp4, HOOp5, HOOp6, HOOp7];
labelStitched = ["Stitched 111","Stitched 110","Stitched 011","Stitched 101","Stitched 100","Stitched 010","Stitched 001"];
imgStitched   = Array{Complex{T},3}(undef, nX, nY, length(HOOps));
weights        = [weightStitched, weightStitched, weightStitched, weightNominal, weightNominal, weightStitched, weightNominal];
idxs           = collect(1:length(HOOps));
for (idx, HOOp, label, weight) in zip(idxs, HOOps, labelStitched, weights)
    @info "[$(idx)] $(label)"
    @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
    imgStitched[:, :, idx] = x;
    fig = plt_image(abs.(x); title=label, vmaxp=99.9)
end


# 1. standard "111"
HOOp1 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("111"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 2. standard "110"
HOOp2 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("110"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 3. standard "011"
HOOp3 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("011"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 4. standard "101"
HOOp4 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("101"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 5. standard "100"
HOOp5 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("100"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 6. standard "010"
HOOp6 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("010"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# 7. standard "001"
HOOp7 = HighOrderOp(gridding, T.(ksphaStandard[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("001"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);

HOOps          = [HOOp1, HOOp2, HOOp3, HOOp4, HOOp5, HOOp6, HOOp7];
labelStandard  = ["Standard 111", "Standard 110", "Standard 011", "Standard 101", "Standard 100", "Standard 010", "Standard 001"];
imgStandard    = Array{Complex{T},3}(undef, nX, nY, length(HOOps));
weights        = [weightStandard, weightStandard, weightStandard, weightNominal, weightNominal, weightStandard, weightNominal];
idxs           = collect(1:length(HOOps));
for (idx, HOOp, label, weight) in zip(idxs, HOOps, labelStandard, weights)
    @info "[$(idx)] $(label)"
    @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
    imgStandard[:, :, idx] = x;
    fig = plt_image(abs.(x); title=label, vmaxp=99.9)
end

# 1. Nominal
labelNominal = ["Nominal 010",];
imgNominal   = Array{Complex{T},3}(undef, nX, nY, 1);

HOOp = HighOrderOp(gridding, T.(ksphaNominal[:, 1:9]'), T.(datatime); sim_method=BlochHighOrder("010"), tr_nominal=ksphaNominal[:, 2:4]', 
        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weightNominal), recParams);
plt_image(abs.(x), vmaxp=99.9)
imgNominal[:, :, 1] = x;

recon = Dict("imgStitched" => imgStitched, "imgStandard" => imgStandard, "imgNominal" => imgNominal,
             "labelStitched" => labelStitched, "labelStandard" => labelStandard, "labelNominal" => labelNominal,);

data["recon"] = recon;
MAT.matwrite(data_file, data)
