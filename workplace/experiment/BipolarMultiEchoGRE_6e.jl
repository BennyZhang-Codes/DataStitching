using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot
import KomaHighOrder.MRIBase: AcquisitionData

path       = "/home/jyzhang/Desktop/pulseq/20240902_invivo/"
seq_file   = path * "grec6e_fov200_200_bw556.seq"
raw_file   = path * "meas_MID00039_FID47484_pulseq_v0_grec6e.mrd"
T          = Float32

seq          = read_seq(seq_file);
TE = seq.DEF["TE"];  # s
ReadoutMode = seq.DEF["ReadoutMode"]; # "Bipolar" or "Monopolar"
# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(raw_file))

shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape)
kdata = get_kdata(raw, shape);
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];
println(size(kdata)); println(kdims);

if ReadoutMode=="Bipolar"
    kdata[:,:,:, collect(range(2, nCon, step=2))] = kdata[:,:,end:-1:1, collect(range(2, nCon, step=2))]; # reverse the even echoes because of bipolar readout
end


############
# IFFT
############
imgs = convert_ifft(kdata, dims=[2,3]);
imgs = CoilCombineSOS(abs.(imgs), 1);
imgs = permutedims(imgs, [3,1,2]);
fig = plt_images(imgs,width=7.5, height=5, vminp=0, vmaxp=99)
fig.savefig("$(path)/$(raw.params["protocolName"])_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)



############
# espirit
############
# using firtst echo for coil sensitivity estimation
Con = 1

_, ktraj_adc = get_kspace(seq);

ktraj = reshape(ktraj_adc[:,1:2], nX, nCon*nY, :)[:, collect(range(Con,nCon*nY,step=nCon)), :];
ktraj = reshape(ktraj.*0, nX*nY, :);
tr    = Trajectory(T.(ktraj'), nY, nX, circular=false, cartesian=true); #plot_traj2d(tr)

k_img = Array{Complex{T},2}(undef,nX*nY,nCha);
for y = 1:nY
    k_img[(y-1)*nX+1:y*nX,:] = transpose(kdata[:, y, : ,Con])
end
# plt_image(abs.(kdata[1, :, : ,Con]).^0.02)
# fig, ax = plt.subplots(1,1); ax.plot(abs.(k_img[:,1]));

dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(k_img,:,nCha);
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));

# espirit
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)

smap = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nY, nX, nCha, 1)

fig = plt_images(permutedims(abs.(smap), [3,1,2]),width=10, height=5)
fig.savefig("$(path)/$(raw.params["protocolName"])_Con$(Con)_CoilSens_espirit.png", dpi=300, bbox_inches="tight", pad_inches=0)


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
images = convert_ifft(kdata[:,:,:,:], dims=[2,3]);   # (nCha, nY, nX, nCon)
ydata = permutedims(images, [2,3,1,4]);    # (nY, nX, nCha, nCon)

# MRIFieldmaps.jl - 3. NCG: diagonal preconditioner

yik_sos = sum(conj(smap) .* ydata; dims=3); # coil combine
yik_sos = yik_sos[:,:,1,:]; # (nY, nX, nCon)
fig = plt_images(permutedims(abs.(yik_sos), [3,1,2]),width=7.5, height=5)
fig = plt_images(permutedims(angle.(yik_sos), [3,1,2]),width=7.5, height=5)
fig.savefig("$(path)/$(raw.params["protocolName"])_pha.png", dpi=300, bbox_inches="tight", pad_inches=0)

(yik_sos_scaled, scale) = b0scale(yik_sos, echotime);
# fig = plt_images(permutedims(abs.(yik_sos_scaled), [3,1,2]),width=10, height=5)

# b0 init
finit = b0init(ydata, echotime; smap); plt_image(finit)

yik_scale = ydata / scale;
fmap_run = (niter, precon, track; kwargs...) -> b0map(yik_scale, echotime; smap, mask,
                            order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

function runner(niter, precon; kwargs...)
    (fmap, times, out) = fmap_run(niter, precon, true; kwargs...) # tracking run
    return (fmap, out.fhats, out.costs, times)
end;

niter_cg_d = 100
(fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag);
b0 = fmap_cg_d * -1;  # Hz


fig = plt_B0map(b0, width=5, height=4)
fig.savefig("$(path)/$(raw.params["protocolName"])_b0map.png", dpi=300, bbox_inches="tight", pad_inches=0, transparent=true)