using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData
using Statistics

T = Float32
#=
data format:
{resolution}_{matrixsize}_{undersamplingrate}, e.g., "0p5_400_r4"

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
out_path     = "$(path)/gre"
if ispath(out_path) == false mkpath(out_path) end

gre_seqs  = ["gres6e_fov200_200_tr25_fa10_bw556.seq",
             "gres6e_fov200_280_tr25_fa10_bw595.seq",
             "gres6e_fov200_332_tr25_fa10_bw602.seq",
             "gres6e_fov200_400_tr25_fa10_bw833.seq"]
gre_mrds  = ["meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mrd",
             "meas_MID00122_FID53010_pulseq_v0_gres6_0p71_standard.mrd",
             "meas_MID00121_FID53009_pulseq_v0_gres6_0p6_standard.mrd",
             "meas_MID00128_FID53016_pulseq_v0_gres6_0p5_standard.mrd"]
filenames = ["gre_1p0_200_r1",
             "gre_0p71_280_r1",
             "gre_0p6_332_r1",
             "gre_0p5_400_r1"]

for (gre_seq, gre_mrd, filename) in zip(gre_seqs, gre_mrds, filenames)

    gre_seq_file = "$(path)/seq/$(gre_seq)" 
    gre_mrd_file = "$(path)/mrd/$(gre_mrd)"


    kdata, ktraj, kdims, shape, TE, ReadoutMode, acqData, raw = read_gre(gre_seq_file, gre_mrd_file);
    nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;

    # # IFFT of gre
    # gre_imgs = convert_ifft(kdata, dims=[2,3]);
    # gre_imgs = CoilCombineSOS(abs.(gre_imgs), 1);
    # gre_imgs = permutedims(gre_imgs, [3,1,2]);
    # fig = plt_images(gre_imgs,width=5, height=5, vminp=0, vmaxp=99)
    # fig.savefig("$(out_path)/$(basename(gre_mrd_file)[1:end-4])_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)
    # MAT.matwrite("$(out_path)/gres6e_1p0.mat", Dict("gre_imgs" => gre_imgs))

    ############
    # espirit
    ############
    sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.987);  # (nX, nY, 1, nCha)
    csm = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nX, nY, 1, nCha) => (nY, nX, nCha, 1) => (nY, nX, nCha)
    fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=5, height=5)  # (nCha, nY, nX)
    fig.savefig("$(out_path)/$(filename)_csm.png", dpi=300, bbox_inches="tight", pad_inches=0)

    ############
    # ΔB₀
    ############
    TE = TE;  # s
    smap = csm[:,:,:]; # [1,2,3,5,6,7,8] (nY, nX, nCha)
    # create mask from Coil-Sensitivity Map
    mask = get_mask(smap[:,:,1], threshold=0); 
    fig = plt_image(mask; vmax=1, vmin=0, width=5, height=5)
    fig.savefig("$(out_path)/$(filename)_mask.png"  , dpi=300, bbox_inches="tight", pad_inches=0)

    images = convert_ifft(kdata[:,:,:,:], dims=[2,3]); # (nCha, nY, nX, nCon)
    ydata = permutedims(images, [2,3,1,4]);            # (nCha, nY, nX, nCon) => (nY, nX, nCha, nCon)

    b0, yik_sos = get_B0map(ydata, TE, smap, mask);
    fig = plt_images(  abs.(yik_sos); dim=3, width=5, height=5, vminp=0, vmaxp=99.0)
    fig.savefig("$(out_path)/$(filename)_mag.png"  , dpi=300, bbox_inches="tight", pad_inches=0)

    fig = plt_images(angle.(yik_sos); dim=3, width=5, height=5, vmin=-pi, vmax=pi)
    fig.savefig("$(out_path)/$(filename)_pha.png"  , dpi=300, bbox_inches="tight", pad_inches=0)

    fig = plt_B0map(b0; width=6, height=5, vmin=-150, vmax=150, color_facecolor="#1f1f1f")
    fig.savefig("$(out_path)/$(filename)_db0.png", dpi=300, bbox_inches="tight", pad_inches=0)

    MAT.matwrite("$(out_path)/$(filename).mat", Dict("img" => yik_sos, "csm" => csm, "b0" => b0, "mask" => mask))
end
