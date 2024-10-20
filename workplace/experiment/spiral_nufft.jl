using KomaHighOrder, MRIReco, PlotlyJS
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory, subsampleIndices, rawdata
import KomaHighOrder.MRIBase: uniqueidx, encSteps1, encSteps2


path       = "/home/jyzhang/Desktop/pulseq/20240828/phantom/"
seq_file   = path * "xw_sp2d_7T-1mm-r4-4segs.seq"
raw_file   = path * "meas_MID00169_FID47132_pulseq_v0_spiral_r4.mrd"
name       = "xw_sp2d_7T-1mm-r4-4segs"
T          = Float64

# get trajectory data
seq          = read_seq(seq_file);  # plot_seq(seq)
_, ktraj_adc = get_kspace(seq[60:end-3]);

Nx, Ny, Nz = Int.(seq.DEF["matrixSize"])
Rep = 1

# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(raw_file))
raw.params["trajectory"] = "custom";

shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape
println(shape)
kdata = Complex{T}.(get_kdata(raw, shape));
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [dims[idx] for idx in 1:length(shape) if shape[idx]>1];


tr           = Trajectory(T.(ktraj_adc[:,1:2]'), nSet, nX, circular=false);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(transpose(reshape(kdata[:,:,Rep,:], 32, nX*nSet)),:,nCha);
acqData      = AcquisitionData(tr, dat, encodingSize=(Nx, Ny));

# acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C

recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
rec      = reconstruction(acqData, recParams);
images   = reshape(rec.data, Nx, Ny, :);
images   = (abs.(images) * prod(size(images)[1:2]));

image_mag = CoilCombineSOS(images, 3)
# plt_image(rotr90(image_mag))
fig = plt_image(image_mag, width=5, height=5)
fig.tight_layout(pad=0, h_pad=0, w_pad=0)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)

# all  channels
fig = plt_images(permutedims(images, [3,1,2]),width=10, height=5)
fig.tight_layout(pad=0, h_pad=0, w_pad=0)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_32cha.png", dpi=300, bbox_inches="tight", pad_inches=0)
