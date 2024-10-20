using KomaHighOrder, MRIReco, PlotlyJS
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory, subsampleIndices, rawdata
import KomaHighOrder.MRIBase: uniqueidx, encSteps1, encSteps2


path       = "/home/jyzhang/Desktop/pulseq/20240924_skope_exp/dat_phantom_big/"
seq_file   = "$(path)/seq/xw_sp2d_7T-1mm-200-r4-noSync-fa10.seq"
mrd_file   = "$(path)/mrd/meas_MID00062_FID49266_pulseq_v0_seg4_r4.mrd"

# path       = "/home/jyzhang/Desktop/pulseq/20240902_invivo/"
# seq_file   = "$(path)/xw_sp2d_7T-1mm-200-r4-4segs.seq"
# mrd_file   = "$(path)/meas_MID00038_FID47483_pulseq_v0_r4.mrd"

T          = Float64

# get trajectory data
seq         = read_seq(seq_file)[end-9:end-3];  # plot_seq(seq)
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)

_, ktraj_adc = get_kspace(seq);

Nx, Ny, Nz = Int.(seq.DEF["matrixSize"])
Rep = 4

# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(mrd_file))
raw.params["trajectory"] = "custom";

shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape
println(shape)
kdata = Complex{T}.(get_kdata(raw, shape));
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];


tr           = Trajectory(T.(ktraj_adc[:,1:2]'), nSet, nX, circular=false);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
# dat[1,1,1]   = reshape(transpose(reshape(kdata[:,:,Rep,:], 32, nX*nSet)),:,nCha);
dat[1,1,1]   = permutedims(reshape(kdata, nCha, nX*nSet), [2,1]);

acqData      = AcquisitionData(tr, dat, encodingSize=(Nx, Ny));

# acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]));  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C;

recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
rec      = reconstruction(acqData, recParams);
images   = reshape(rec.data, Nx, Ny, :);
images   = (abs.(images) * prod(size(images)[1:2]));

image_mag = CoilCombineSOS(images, 3);
fig = plt_image(rotr90(image_mag, 1), width=5, height=5)
fig.tight_layout(pad=0, h_pad=0, w_pad=0)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_sos111.png", dpi=300, bbox_inches="tight", pad_inches=0)

# all  channels
imgs = zeros(Ny, Nx, nCha);
for i in 1:nCha
    imgs[:,:,i] = rotr90(images[:,:,i]);
end
fig = plt_images(permutedims(imgs, [3,1,2]),width=10, height=5)
fig.tight_layout(pad=0, h_pad=0, w_pad=0)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_32cha.png", dpi=300, bbox_inches="tight", pad_inches=0)


solver = "cgnr"
reg = "L2"
iter = 2
λ = 1e-3
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity[end:-1:1,:,:,:], Nx, Ny, 1, nCha)); # reverse the x-axis

# recParams[:correctionMap] = 1im.*Matrix(transpose(b0)) .* 2π;  

@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)



# HighOrderOp
_, ktraj_adc = get_kspace(seq);
times = KomaMRIBase.get_adc_sampling_times(seq);
tr_nominal = Trajectory(ktraj_adc'[1:3,:], nSet, nX; circular=false, times=times);
tr_dfc     = Trajectory(zeros(9, nX*nSet), nSet, nX; circular=false, times=times);
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc , BlochHighOrder("000"); fieldmap=Matrix(rotl90(b0)), Nblocks=9, grid=1);

smaps = sensitivity[end:-1:1,:,1,:];
S = SensitivityOp(reshape(Complex{T}.(smaps),:,nCha),1);
Op = DiagOp(Op, nCha) ∘ S ;
# recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCha));
recParams[:encodingOps] = reshape([Op], 1,1);

@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

