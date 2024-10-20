using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot
import KomaHighOrder.MRIBase: AcquisitionData

path       = "/home/jyzhang/Desktop/pulseq/20240924_skope_exp/dat_phantom_acr_isocenter/"
seq_file   = "$(path)/seq/gres6e_fov200_200_bw556.seq"
raw_file   = "$(path)/mrd/meas_MID00071_FID49275_pulseq_v0_gres6.mrd"
T          = Float64

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
fig.savefig("$(path)/$(raw.params["protocolName"])_Con1_CoilSens_espirit.png", dpi=300, bbox_inches="tight", pad_inches=0)

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




# path = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo/"
seq_file = "$(path)/seq/xw_sp2d_7T-1mm-200-r4-noSync-fa10.seq"
mrd_file = "$(path)/mrd/meas_MID00074_FID49278_pulseq_v0_seg1_r4.mrd"
dfc_file = "$(path)/dfc/xw_sp2d_7T-1mm-200-r4.mat"


########################################################################
# 1. load the *.seq file and extract some infomations
########################################################################
seq = read_seq(seq_file)[end-9:end-3]; 
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
_, ktrNominal = get_kspace(seq);
times = KomaMRIBase.get_adc_sampling_times(seq);
Nx, Ny, Nz = Int.(seq.DEF["matrixSize"])


########################################################################
# 2. load the *.mrd file and extract the k-space signal 
########################################################################
raw = RawAcquisitionData(ISMRMRDFile(mrd_file))
raw.params["trajectory"] = "custom";
shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape);
kdata = Complex{T}.(get_kdata(raw, shape));
kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];

tr           = Trajectory(T.(ktrNominal[:,1:2]'), nSet, nX, circular=false);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = permutedims(reshape(kdata, nCha, nX*nSet), [2,1]);
acqData      = AcquisitionData(tr, dat, encodingSize=(Nx, Ny));
C = maximum(2*abs.(acqData.traj[1].nodes[:]));  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C;

########################################################################
# 3. load the *.mat file and extract the DFC data
########################################################################
dfcStitched, dfcStandard, ntStitched, ntStandard = load_dfc_mat(dfc_file);
hoseqStitched = HO_Sequence(seq);                     # hoseq, defined in KomaHighOrder.jl
hoseqStandard = HO_Sequence(seq);                     # hoseq, defined in KomaHighOrder.jl
hoseqStitched.GR_dfc[2:4, :] = hoseqStitched.SEQ.GR;  # copy the 1st-order gradient data from the seq object to the hoseq object
hoseqStandard.GR_dfc[2:4, :] = hoseqStandard.SEQ.GR;  # copy the 1st-order gradient data from the seq object to the hoseq object
hoseqStitched.GR_dfc[:,6] = dfcStitched;              # "6" is the index of the readout block in the spiral sequence
hoseqStandard.GR_dfc[:,6] = dfcStandard;              # "6" is the index of the readout block in the spiral sequence
plot_seq(hoseqStitched)
plot_seq(hoseqStandard)
# finally, hoseq* contains both the nominal trajectory and the measured trajectory (up to 2nd-order)


########################################################################
# 4. load the *.mat file and extract the DFC data
########################################################################
_, K_nominal, _, K_dfcStitched = get_kspace(hoseqStitched; Δt=1);
_, _, _, K_dfcStandard = get_kspace(hoseqStandard; Δt=1);


tr_nominal      = Trajectory(  K_nominal'[1:3,:], nSet, nX; circular=false, times=times);
tr_dfc_stitched = Trajectory(K_dfcStitched'[:,:], nSet, nX; circular=false, times=times);
tr_dfc_standard = Trajectory(K_dfcStandard'[:,:], nSet, nX; circular=false, times=times);



########################################################################
# NUFFT
########################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
rec      = reconstruction(acqData, recParams);
images   = reshape(rec.data, Nx, Ny, :);
# images   = (abs.(images) * prod(size(images)[1:2]));
# all  channels
imgs = Array{ComplexF32,3}(undef, Ny, Nx, nCha);
for i in 1:nCha
    imgs[:,:,i] = rotr90(images[:,:,i]);
end

img = sum(conj(smap) .* imgs; dims=3)[:,:,1]; # plt_image(abs.(img))
# image_mag = CoilCombineSOS(images, 3);
fig = plt_image(abs.(img), width=5, height=5)
fig.tight_layout(pad=0, h_pad=0, w_pad=0)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_sos.png", dpi=300, bbox_inches="tight", pad_inches=0)


fig = plt_images(permutedims(abs.(imgs), [3,1,2]),width=10, height=5)
fig.tight_layout(pad=0, h_pad=0, w_pad=0)
fig.savefig("$(path)/$(raw.params["protocolName"])_nufft_32cha.png", dpi=300, bbox_inches="tight", pad_inches=0)


########################################################################
# SENSE
########################################################################
solver = "cgnr"
reg = "L2"
iter = 30
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

########################################################################
# HighOrderOp
########################################################################
_, ktraj_adc = get_kspace(seq);
times = KomaMRIBase.get_adc_sampling_times(seq);
tr_nominal = Trajectory(ktraj_adc'[1:3,:], nSet, nX; circular=false, times=times);
tr_dfc     = Trajectory(zeros(9, nX*nSet), nSet, nX; circular=false, times=times);
Op0 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc , BlochHighOrder("000"); fieldmap=Matrix(rotl90(b0)), Nblocks=9, grid=1);

B0map = rotl90(b0); Nblocks=2;

solver = "cgnr"
reg = "L2"
iter = 10
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

smaps = sensitivity[end:-1:1,:,1,:];
S = SensitivityOp(reshape(Complex{T}.(smaps),:,nCha),1);

# 1. nominal trajectory, BlochHighOrder("000")
Op0 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1, verbose=false);
Op = DiagOp(Op0, nCha) ∘ S ;
recParams[:encodingOps] = reshape([Op], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_nominal_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 2. stitched trajectory, BlochHighOrder("111")
Op0 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1, verbose=false);
Op = DiagOp(Op0, nCha) ∘ S ;
recParams[:encodingOps] = reshape([Op], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_stitched111_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 3. stitched trajectory, BlochHighOrder("110")
Op0 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1, verbose=false);
Op = DiagOp(Op0, nCha) ∘ S ;
recParams[:encodingOps] = reshape([Op], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_stitched110_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 4. standard trajectory, BlochHighOrder("111")
Op0 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1, verbose=false);
Op = DiagOp(Op0, nCha) ∘ S ;
recParams[:encodingOps] = reshape([Op], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_standard111_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 5. standard trajectory, BlochHighOrder("110")
Op0 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_standard , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=Matrix(B0map), grid=1, verbose=false);
Op = DiagOp(Op0, nCha) ∘ S ;
recParams[:encodingOps] = reshape([Op], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotr90(rec))
fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_standard110_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

