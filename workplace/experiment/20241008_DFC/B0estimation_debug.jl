# https://github.com/HarmonizedMRI/Calibration/blob/main/b0/b0estimation/b0estimation.jl
using MRIFieldmaps: b0init, b0map
using ROMEO: unwrap
using MAT: matread, matwrite
using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot
import KomaHighOrder.MRIBase: AcquisitionData

path       = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo"
seq_file = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("gres6", f)][1]
raw_file = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin("gres6", f)][1]
T          = Float64

kdata, ktraj, kdims, shape, TE, ReadoutMode, gre_acqData, raw = read_gre(seq_file, raw_file);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;

# IFFT of gre
gre_imgs = convert_ifft(kdata, dims=[2,3]);
gre_imgs = CoilCombineSOS(abs.(gre_imgs), 1);
gre_imgs = permutedims(gre_imgs, [3,1,2]);
fig = plt_images(gre_imgs,width=7.5, height=5, vminp=0, vmaxp=99)

############
# espirit
############
sensitivity = espirit(gre_acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)
csm = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nX, nY, 1, nCha) => (nY, nX, nCha, 1) => (nY, nX, nCha)
fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=10, height=5)  # (nCha, nY, nX)

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
fig = plt_B0map(b0, width=5, height=4)



# MRIFieldmaps.jl - 3. NCG: diagonal preconditioner
yik_sos = sum(conj(smap) .* ydata; dims=3)[:,:,1,:]; # coil combine => (nY, nX, 1, nCon) => (nY, nX, nCon)

(yik_sos_scaled, scale) = b0scale(yik_sos, TE); # fig = plt_images(permutedims(abs.(yik_sos_scaled), [3,1,2]),width=10, height=5)

# b0 init
# finit = b0init(ydata, TE; smap); plt_image(finit; title="b0init")

yik_scale = ydata / scale;
fmap_run = (niter, precon, track; kwargs...) -> b0map(yik_scale, TE; smap, mask, order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

function runner(niter, precon; kwargs...)
    (fmap, times, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
    return (fmap, out.fhats, out.costs, times)
end;

niter_cg_d =200
(fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag);
b0 = fmap_cg_d * -1;  # Hz
plt_B0map(b0, width=5, height=4)