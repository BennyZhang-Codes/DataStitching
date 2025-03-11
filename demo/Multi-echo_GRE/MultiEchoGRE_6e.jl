using HighOrderMRI, MRIReco, MRICoilSensitivities
import HighOrderMRI.MRIBase: AcquisitionData

path     = joinpath(@__DIR__, "demo/Multi-echo_GRE")
seq_file = "$(path)/gres6e_1p0_200_r1_tr25_fa10_bw556.seq"               # *.seq file is the pulseq's sequence file
mrd_file = "$(path)/meas_MID00117_FID53005_pulseq_gres6_1p0_200_r1.mrd"  # *.mat file contains the dynamic field data from both stitching method and the standard method.


seq         = read_seq(seq_file);
TE          = seq.DEF["TE"];  # s
ReadoutMode = seq.DEF["ReadoutMode"]; # "Bipolar" or "Monopolar" or "Separate"
raw         = RawAcquisitionData(ISMRMRDFile(mrd_file)) # get signal data

shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape)
kdata = get_kdata(raw, shape);
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];

if ReadoutMode=="Bipolar"
    kdata[:,:,:, collect(range(2, nCon, step=2))] = kdata[:,:,end:-1:1, collect(range(2, nCon, step=2))]; # reverse the even echoes because of bipolar readout
end

@info TE, ReadoutMode
@info size(kdata), kdims

############
# IFFT
############
imgs = convert_ifft(kdata, dims=[2,3]);
imgs = CoilCombineSOS(abs.(imgs), 1);
imgs = permutedims(imgs, [3,1,2]);
fig  = plt_images(imgs; width=5, height=5, vminp=0, vmaxp=99)

############
# espirit
############
# using firtst echo for coil sensitivity estimation
Con = 1
T   = Float32
_, ktraj_adc = get_kspace(seq);
ktraj = reshape(ktraj_adc[:,1:2], nX, nCon*nY, :)[:, collect(range(Con,nCon*nY,step=nCon)), :];
ktraj = reshape(ktraj.*0, nX*nY, :);
tr    = Trajectory(T.(ktraj'), nY, nX, circular=false, cartesian=true); #plot_traj2d(tr)
k_img = Array{Complex{T},2}(undef,nX*nY,nCha);
for y = 1:nY
    k_img[(y-1)*nX+1:y*nX,:] = transpose(kdata[:, y, : ,Con])
end

dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(k_img,:,nCha);
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));

# espirit
sensitivity  = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)
smap         = permutedims(sensitivity, [2,1,4,3])[:,:,:,1]; # (nY, nX, nCha, 1)
fig          = plt_images(permutedims(abs.(smap), [3,1,2]); width=5, height=5)


############
# ΔB₀
############
using MRIFieldmaps: b0scale, b0init, b0map
# mask from Coil-Sensitivity Map
mask = smap[:,:,1];
mask[abs.(mask) .> 0] .= 1;
mask = isone.(mask); 
plt_image(mask)

echotime = seq.DEF["TE"];  # s
images   = convert_ifft(kdata[:,:,:,:], dims=[2,3]);   # (nCha, nY, nX, nCon)
ydata    = permutedims(images, [2,3,1,4]);    # (nY, nX, nCha, nCon)

# MRIFieldmaps.jl - 3. NCG: diagonal preconditioner

yik_sos = sum(conj(smap) .* ydata; dims=3); # coil combine
yik_sos = yik_sos[:,:,1,:]; # (nY, nX, nCon)
fig = plt_images(permutedims(  abs.(yik_sos), [3,1,2]); width=5, height=5, vminp=0, vmaxp=99)
fig = plt_images(permutedims(angle.(yik_sos), [3,1,2]); width=5, height=5, vmin=-π, vmax=π)

(yik_sos_scaled, scale_b0) = b0scale(yik_sos, echotime);

# b0 init
finit = b0init(ydata, echotime; smap); 
plt_image(finit)

yik_scale = ydata / scale_b0;
fmap_run  = (niter, precon, track; kwargs...) -> b0map(yik_scale, echotime; smap, mask,
                            order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

function runner(niter, precon; kwargs...)
    (fmap, times, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
    return (fmap, out.fhats, out.costs, times)
end

niter_cg_d = 40
(fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag);
b0 = fmap_cg_d;  # Hz

fig = plt_B0map(b0; width=5, height=5, vmin=-150, vmax=150, cmap="jet")
