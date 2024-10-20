path       = "/home/jyzhang/Desktop/pulseq/20240827/"
seq_gre    = path * "gre6e_fov150_150.seq"
raw_gre    = path * "meas_MID00062_FID46740_pulseq_v0_gre6e.mrd"
T          = Float32
Con = 3

# get trajectory data
seq          = read_seq(seq_gre);
_, ktraj_adc = get_kspace(seq);

# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(raw_gre))
raw.params["trajectory"] = "custom";

shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape
println(shape)
kdata = get_kdata(raw, shape);
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [dims[idx] for idx in 1:length(shape) if shape[idx]>1];
println(size(kdata))
println(kdims)

ktraj = reshape(ktraj_adc[:,1:2], nX, nCon*nY, :)[:, collect(range(Con,nCon*nY,step=nCon)), :];
ktraj = reshape(ktraj, nX*nY, :);
tr           = Trajectory(T.(ktraj'), nY, nX, circular=false, cartesian=true); plot_traj2d(tr);

k_img = Array{Complex{T},2}(undef,nX*nY,nCha);
for y = 1:nY
    k_img[(y-1)*nX+1:y*nX,:] = transpose(kdata[:, y, : ,Con])
end
plt_image(abs.(kdata[1, :, : ,Con]).^0.02)
fig, ax = plt.subplots(1,1)
ax.plot(abs.(k_img[:,1]))


dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(k_img,:,nCha);
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));


# espirit
# acqDataCart = regrid(acqData, (Nx, Ny); cgnr_iter=3);
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98);
plot_imgs_subplots(real(sensitivity[:,:,1,:]), 4, 8)
plot_img(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]))
p_smap_espirit = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), 4, 8; title="Coil Sensitivity (espirit)")
savefig(p_smap_espirit,  path * "$(raw.params["protocolName"])-CoilSens_espirit.svg", format="svg", height=400, width=800)
