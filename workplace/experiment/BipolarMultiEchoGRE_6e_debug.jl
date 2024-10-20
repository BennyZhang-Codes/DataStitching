using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot
import KomaHighOrder.MRIBase: AcquisitionData

path       = "/home/jyzhang/Desktop/pulseq/20240828/invivo/"
seq_file   = path * "gre6e_fov150_150_bw833.seq"
raw_file   = path * "meas_MID00131_FID47094_pulseq_v0_gre6e_833.mrd"
T          = Float32

seq          = read_seq(seq_file);
# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(raw_file))

shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape)
kdata = get_kdata(raw, shape);
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];
println(size(kdata)); println(kdims);

kdata[:,:,:, collect(range(2, nCon, step=2))] = kdata[:,:,end:-1:1, collect(range(2, nCon, step=2))] # reverse the even echoes because of bipolar readout



############
# NUFFT
############
imgs = convert_ifft(kdata, dims=[2,3]);
imgs = CoilCombineSOS(abs.(imgs), 1);
imgs = permutedims(imgs, [3,1,2]);
fig = plt_images(imgs,width=7.5, height=5)
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
plt_image(abs.(kdata[1, :, : ,Con]).^0.02)
fig, ax = plt.subplots(1,1); ax.plot(abs.(k_img[:,1]));

dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(k_img,:,nCha);
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));

# espirit
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)
sensitivity = permutedims(sensitivity, [2,1,4,3]);# (nY, nX, nCha, nCon)

fig = plt_images(permutedims(abs.(sensitivity)[:,:,:,1], [3,1,2]),width=10, height=5)
fig.savefig("$(path)/$(raw.params["protocolName"])_Con$(Con)_CoilSens_espirit.png", dpi=300, bbox_inches="tight", pad_inches=0)


############
# ΔB₀
############
TE = seq.DEF["TE"]
images = convert_ifft(kdata, dims=[2,3]);   # (nCha, nY, nX, nCon)
images = permutedims(images, [2,3,1,4]);    # (nY, nX, nCha, nCon)

fig = plt_images(permutedims(abs.(images)[:,:,:,1], [3,1,2]),width=10, height=5)
fig = plt_images(permutedims(angle.(images)[:,:,:,1], [3,1,2]),width=10, height=5)

yik_sos = sum(conj(sensitivity) .* images; dims=3) # coil combine

fig = plt_images(permutedims(abs.(yik_sos)[:,:,1,:], [3,1,2]),width=10, height=5)
fig = plt_images(permutedims(angle.(yik_sos)[:,:,1,:], [3,1,2]),width=10, height=5)


# estimate_b0(images, TE; sensitivity[:,:,:,1])

data = images
echotimes = TE
smap = sensitivity[:,:,:,1]

finit = b0init(data, echotimes; smap); plt_image(finit)

# Compute the square root sum-of-squares coil-combined image
# to aid with the unwrapping procedure
imsos = sqrt.(dropdims(sum(abs2, data; dims = 3); dims = 3))

# Convert the initial field map to radians
ΔTE = echotimes[2] - echotimes[1]
phase = (2π * ΔTE) .* finit

# Unwrap, then convert back to units of frequency (cycles/unit time)
unwrapped = unwrap(phase; mag = imsos) ./ (2π * ΔTE)
plt_image(unwrapped)
# Regularize the unwrapped field map
(fhat,) = b0map(data, echotimes; finit = unwrapped, smap, mask)

plt_image(fhat)



mask = sensitivity[:,:,1,1];
mask[abs.(mask) .> 0] .= 1;
mask = isone.(mask); plt_image(mask)

yik_sos = sum(conj(sensitivity) .* images; dims=3) # coil combine
(yik_sos_scaled, scale) = b0scale(yik_sos, TE)
yik_scale = images / scale
finit = b0init(images, TE; smap)


b0map_run = (niter, precon, track; kwargs...) ->b0map(yik_scale, TE; finit, smap, mask,
                            order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

function b0map_runner(niter, precon; kwargs...)
    (fmap, times, out) = b0map_run(niter, precon, true; kwargs...)
    return (fmap, out.fhats, out.costs, times)
end;

niter_cg_d = 100  # (400)
(fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = b0map_runner(niter_cg_d, :diag)

b0 = fmap_cg_d * -1  # Hz
plt_image(b0)


# b0 = imresize(b0,size(sensitivity)[1:2], method=BSpline(Cubic()))   # using ImageTransformations, Interpolations


