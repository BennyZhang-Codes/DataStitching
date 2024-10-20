using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot
import KomaHighOrder.MRIBase: AcquisitionData

path       = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo/"
seq_file   = "$(path)/seq/gres6e_fov200_200_bw556.seq"
raw_file   = "$(path)/mrd/meas_MID00066_FID49933_pulseq_v0_gres6.mrd"
T          = Float32


kdata, ktraj, kdims, shape, TE, ReadoutMode, acqData, raw = read_gre(seq_file, raw_file);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;

# IFFT of gre
imgs = convert_ifft(kdata, dims=[2,3]);
imgs = CoilCombineSOS(abs.(imgs), 1);
imgs = permutedims(imgs, [3,1,2]);
fig = plt_images(imgs,width=7.5, height=5, vminp=0, vmaxp=99)
fig.savefig("$(path)/$(raw.params["protocolName"])_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)


############
# espirit
############
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)

csm = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nX, nY, 1, nCha) => (nY, nX, nCha, 1) => (nX, nY, nCha)
fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=10, height=5)  # (nCha, nY, nX)
fig.savefig("$(path)/$(raw.params["protocolName"])_Con$(Con)_CoilSens_espirit.png", dpi=300, bbox_inches="tight", pad_inches=0)


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

fig = plt_images(permutedims(  abs.(yik_sos), [3,1,2]),width=7.5, height=5)
fig = plt_images(permutedims(angle.(yik_sos), [3,1,2]),width=7.5, height=5)
fig.savefig("$(path)/$(raw.params["protocolName"])_pha.png"  , dpi=300, bbox_inches="tight", pad_inches=0)


fig = plt_B0map(b0, width=5, height=4)
fig.savefig("$(path)/$(raw.params["protocolName"])_b0map.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)

