using MAT


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

# FindDelay, using a model-based approach
Δx, Δy, Δz = T.(FOV ./ matrixSize);
nX, nY, nZ = matrixSize;



gridding  = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)

solver = "cgnr"; reg = "L2"; iter = 10; λ = 0.
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

weight = SampleDensity(pha_spha'[2:3,:], (nX, nY));

HOOp = HighOrderOp(gridding, T.(pha_spha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
plt_image(abs.(x1), vmaxp=99.9)


