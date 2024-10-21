using KomaHighOrder, MRIReco, PlotlyJS
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory, subsampleIndices, rawdata
import KomaHighOrder.MRIBase: uniqueidx, encSteps1, encSteps2


path       = "/home/jyzhang/Desktop/pulseq/20240902_invivo/"
seq_file   = path * "xw_sp2d_7T-1mm-200-r4-4segs.seq"
raw_file   = path * "meas_MID00038_FID47483_pulseq_v0_r4.mrd"
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
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];
println(size(kdata)); println(kdims);

tr           = Trajectory(T.(ktraj_adc[:,1:2]'), nSet, nX, circular=false);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(transpose(reshape(kdata[:,:,Rep,:], 32, nX*nSet)),:,nCha);
acqData      = AcquisitionData(tr, dat, encodingSize=(Nx,Ny));

# acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C

recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
rec      = reconstruction(acqData, recParams);
images   = reshape(rec.data, Nx, Ny, :);
images   = (abs.(images) * prod(size(images)[1:2]));

image_mag = CoilCombineSOS(images, 3);
# plt_image(rotr90(image_mag))
fig = plt_image(transpose(image_mag), width=5, height=5)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)
# all  channels
fig = plt_images(permutedims(images, [3,2, 1]),width=10, height=5)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_32Cha.png", dpi=300, bbox_inches="tight", pad_inches=0)

# single channel
# plot_image(images[:,:,1]; )

# s = scatter(x=ktraj_adc'[1,:], y=ktraj_adc'[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)


solver = "cgnr"
reg = "L2"
iter = 20
λ = 1e-2
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCha));

# recParams[:correctionMap] = 1im.*Matrix(transpose(b0)) .* 2π;

@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(transpose(rec))


fig.savefig("$(path)/$(raw.params["protocolName"])_sense_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)


# HighOrderOp
_, ktraj_adc = get_kspace(seq[60:end-3]);
times = KomaMRIBase.get_adc_sampling_times(seq[60:end-3]);
tr_nominal = Trajectory(ktraj_adc'[1:3,:], nSet, nX; circular=false, times=times);
tr_dfc     = Trajectory(zeros(9, nX*nSet), nSet, nX; circular=false, times=times);
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc , BlochHighOrder("000"); fieldmap=Matrix(transpose(b0)), Nblocks=9, grid=1);

smaps = sensitivity[:,:,1,:];
S = SensitivityOp(reshape(Complex{T}.(smaps),:,nCha),1);
Op = DiagOp(Op, nCha) ∘ S ;
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCha));
recParams[:encodingOps] = reshape([Op], 1,1);

@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(transpose(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)



# L1-Wavelet regularized CS reconstruction
cs_solver = "admm"
cs_reg    = "L1"
cs_sparse = "Wavelet"
cs_λ      = 2.e-1
cs_iter   = 1000


CSParams = Dict{Symbol, Any}()
CSParams[:reco] = "multiCoil"
CSParams[:reconSize] = (Nx, Ny)
CSParams[:solver] = cs_solver
CSParams[:regularization] = cs_reg
CSParams[:sparseTrafo] = cs_sparse
CSParams[:λ] = cs_λ
CSParams[:iterations] = cs_iter
CSParams[:ρ] = 0.1
CSParams[:absTol] = 1.e-15
CSParams[:relTol] = 1.e-3
CSParams[:tolInner] = 1.e-2
CSParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCha));
CSParams[:normalizeReg] = true
@time rec = abs.(reconstruction(acqData, CSParams).data[:,:]);
fig = plt_image(transpose(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_CS_$(cs_solver)_$(cs_reg)_$(cs_sparse)_$(cs_iter)_$(cs_λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# HighOrderOp
_, ktraj_adc = get_kspace(seq[60:end-3]);
times = KomaMRIBase.get_adc_sampling_times(seq[60:end-3]);
tr_nominal = Trajectory(ktraj_adc'[1:3,:], nSet, nX; circular=false, times=times);
tr_dfc     = Trajectory(zeros(9, nX*nSet), nSet, nX; circular=false, times=times);
Op = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc , BlochHighOrder("000"); fieldmap=Matrix(transpose(b0)), Nblocks=9, grid=1);

smaps = sensitivity[:,:,1,:];
S = SensitivityOp(reshape(Complex{T}.(smaps),:,nCha),1);
Op = DiagOp(Op, nCha) ∘ S ;
CSParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, nCha));
CSParams[:encodingOps] = reshape([Op], 1,1);

@time rec = abs.(reconstruction(acqData, CSParams).data[:,:]);
fig = plt_image(transpose(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_CS_HOO_$(cs_solver)_$(cs_reg)_$(cs_sparse)_$(cs_iter)_$(cs_λ).png", dpi=300, bbox_inches="tight", pad_inches=0)


