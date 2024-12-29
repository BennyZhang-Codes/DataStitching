using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData
using Statistics
#=
prepare data for synchronization
1. Multi-Echo GRE (ME-GRE) data: 
    (1) Coil-Sensitivity Map (CSM)
    (2) ΔB₀ map
2. Spiral imaging datas 
    (1) MRI Signal data

Notes: directory structure:
    - path/seq/*.seq     # sequence files
    - path/mrd/*.mrd     # MRD files
    - path/syn/*.mat     # output directory
=#



path         = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/"
gre_seq_file = "$(path)/seq/gres6e_fov200_400_tr25_fa10_bw833.seq" 
gre_mrd_file = "$(path)/mrd/meas_MID00128_FID53016_pulseq_v0_gres6_0p5_standard.mrd"

if ispath("$(path)/syn") == false mkpath("$(path)/syn") end

T = Float32

kdata, ktraj, kdims, shape, TE, ReadoutMode, acqData, raw = read_gre(gre_seq_file, gre_mrd_file);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;

# IFFT of gre
gre_imgs = convert_ifft(kdata, dims=[2,3]);
gre_imgs = CoilCombineSOS(abs.(gre_imgs), 1);
gre_imgs = permutedims(gre_imgs, [3,1,2]);
fig = plt_images(gre_imgs,width=5, height=5, vminp=0, vmaxp=99)
fig.savefig("$(path)/syn/$(basename(gre_mrd_file)[1:end-4])_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)
MAT.matwrite("$(path)/syn/gres6e_1p0.mat", Dict("gre_imgs" => gre_imgs))

############
# espirit
############
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.987);  # (nX, nY, 1, nCha)
csm = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nX, nY, 1, nCha) => (nY, nX, nCha, 1) => (nY, nX, nCha)
fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=5, height=5)  # (nCha, nY, nX)
fig.savefig("$(path)/syn/$(basename(gre_mrd_file)[1:end-4])_CoilSens_espirit.png", dpi=300, bbox_inches="tight", pad_inches=0)

############
# ΔB₀
############
TE = TE;  # s
smap = csm[:,:,:]; # [1,2,3,5,6,7,8] (nY, nX, nCha)
# create mask from Coil-Sensitivity Map
mask = get_mask(smap[:,:,1], threshold=0); plt_image(mask)
images = convert_ifft(kdata[:,:,:,:], dims=[2,3]); # (nCha, nY, nX, nCon)
ydata = permutedims(images, [2,3,1,4]);            # (nCha, nY, nX, nCon) => (nY, nX, nCha, nCon)

b0, yik_sos = get_B0map(ydata, TE, smap, mask);
fig = plt_images(permutedims(  abs.(yik_sos), [3,1,2]),width=5, height=5)
fig = plt_images(permutedims(angle.(yik_sos), [3,1,2]),width=5, height=5)
fig.savefig("$(path)/syn/$(basename(gre_mrd_file)[1:end-4])_pha.png"  , dpi=300, bbox_inches="tight", pad_inches=0)
fig = plt_B0map(-b0, width=5, height=4)
fig.savefig("$(path)/syn/$(basename(gre_mrd_file)[1:end-4])_b0map.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)

MAT.matwrite("$(path)/syn/syn_$(basename(gre_mrd_file)[1:end-4]).mat", Dict("csm" => csm, "b0" => b0, "mask" => mask))




extract_data(seq::String, mrd::String, path::String) = begin
    seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin(seq, f)][1]
    mrd_file     = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin(mrd, f)][1]
    println(seq_file, mrd_file)
    seq = read_seq(seq_file)[end-9:end-3]; 
    # times = KomaMRIBase.get_adc_sampling_times(seq);
    FOV        = seq.DEF["FOV"]
    matrixSize = Int.(seq.DEF["matrixSize"])

    raw = RawAcquisitionData(ISMRMRDFile(mrd_file))
    raw.params["trajectory"] = "custom";
    shape = get_ksize(raw);
    nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape);
    kdata = Complex{T}.(get_kdata(raw, shape));

    # kdata = mean(kdata, dims=9);
    kdata = kdata[:,:,:,:,:,:,:,:,1,:,:];

    kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
    kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];

    data = permutedims(reshape(kdata, nCha, nX*nSet), [2,1]);
    MAT.matwrite("$(path)/syn/syn_$(basename(mrd_file)[1:end-4]).mat", Dict("data" => data, "FOV" => FOV, "matrixSize" => matrixSize))
end

seqs = ["7T_1mm-200-r4_max51-fa90.seq",
        "7T_1mm-200-r3_max51-fa90.seq",
        "7T_1mm-200-r2_max51-fa90.seq",
        "7T_0.71mm-280-r2_max51-fa90.seq",
        "7T_0.6mm-332-r3_max51-fa90.seq",
        "7T_0.5mm-400-r4_max51-fa90.seq"]
mrds = ["r4_1p0_standard",
        "r3_1p0_standard",
        "r2_1p0_standard",
        "r2_0p71_standard",
        "r3_0p6_standard",
        "r4_0p5_standard"]

for idx = 1:6
    extract_data(seqs[idx], mrds[idx], path)
end

