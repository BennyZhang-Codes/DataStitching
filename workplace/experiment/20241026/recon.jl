using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

T = Float64

path         = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo"
r=4
seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("r$(r)", f)][1]
mrd_file     = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin("r$(r)", f)][1]
syn_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin("r$(r)", f)][1]
gre_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(r"^syn.*gre.*mat$", f)][1]

seq = read_seq(seq_file)[end-9:end-3]; 
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
TE = seq.DEF["TE"];
_, k_adc = get_kspace(seq);
times = KomaMRIBase.get_adc_sampling_times(seq);
times = times .- times[1] .+ TE;

# Coil-Sensitivity Map (CSM), ΔB0 map, mask
b0map      = matread(gre_file)["b0"];
csm        = matread(gre_file)["csm"];    # (nY, nX, nCha)
mask       = matread(gre_file)["mask"];

# MRI signal data, FOV, matrix size, synchronized dynamic fields (kspha, rad, rad/m⁻¹, rad/m⁻²)
data       = matread(syn_file)["data"];
FOV        = matread(syn_file)["FOV"];
matrixSize = matread(syn_file)["matrixSize"];
kspha      = matread(syn_file)["ksphaStandard_syn"]/2π;


nY, nX, nCha = size(csm);
csm = permutedims(csm, [2,1,3]); # (nX, nY, nCha)
csm = Complex{T}.(reshape(csm[end:-1:1,:,:,:], nX, nY, 1, nCha)); # reverse the x-axis


nSample, nCha = size(data);
tr           = Trajectory(T.(kspha[:,2:3]'), 1, nSample, circular=false, times=times);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = data;
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));
C = maximum(2*abs.(acqData.traj[1].nodes[:]));  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C;


tr_nominal = Trajectory(k_adc'[1:3,:], 1, nSample; circular=false, times=times);
tr_kspha   = Trajectory(kspha'[:,:],   1, nSample; circular=false, times=times);

########################################################################
# SENSE
########################################################################
solver = "cgnr"
reg = "L2"
iter = 20
λ = 1e-3
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (nX, nY)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = csm;

@time rec = reconstruction(acqData, recParams).data[:,:];
fig = plt_image(rotr90(abs.(rec)))

fig.savefig("$(path)/$(raw.params["protocolName"])_sense_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)



########################################################################
# HighOrderOp
########################################################################
b0map = Matrix(rotl90(b0map)); 
Nblocks=5;

Op0 = HighOrderOp((nX, nY), tr_nominal, tr_kspha, BlochHighOrder("000"); fieldmap=b0map, Nblocks=9, grid=1);

solver = "cgnr"
reg = "L2"
iter = 10
λ = 1e-3
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (nX, nY)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
# recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = csm;

S = SensitivityOp(reshape(Complex{T}.(csm[:,:,1,:]),:,nCha),1);

# 1. nominal trajectory, BlochHighOrder("000")
Op0 = HighOrderOp((nX, nY), tr_nominal, tr_kspha , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=b0map, grid=1, verbose=false);
recParams[:encodingOps] = reshape([DiagOp(Op0, nCha) ∘ S], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]); fig = plt_image(rotr90(rec))

fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_nominal_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 2. stitched trajectory, BlochHighOrder("111")
Op0 = HighOrderOp((nX, nY), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=b0map, grid=1, verbose=false);
recParams[:encodingOps] = reshape([DiagOp(Op0, nCha) ∘ S], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]); fig = plt_image(rotr90(rec))

fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_stitched111_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 3. stitched trajectory, BlochHighOrder("110")
Op0 = HighOrderOp((nX, nY), tr_nominal, tr_dfc_stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=b0map, grid=1, verbose=false);
recParams[:encodingOps] = reshape([DiagOp(Op0, nCha) ∘ S], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]); fig = plt_image(rotr90(rec))

fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_stitched110_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 4. standard trajectory, BlochHighOrder("111")
Op0 = HighOrderOp((nX, nY), tr_nominal, tr_dfc_standard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=b0map, grid=1, verbose=false);
recParams[:encodingOps] = reshape([DiagOp(Op0, nCha) ∘ S], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]); fig = plt_image(rotr90(rec))

fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_standard111_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)

# 5. standard trajectory, BlochHighOrder("110")
Op0 = HighOrderOp((nX, nY), tr_nominal, tr_dfc_standard , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=b0map, grid=1, verbose=false);
recParams[:encodingOps] = reshape([DiagOp(Op0, nCha) ∘ S], 1,1);
@time rec = abs.(reconstruction(acqData, recParams).data[:,:]); fig = plt_image(rotr90(rec))

fig.savefig("$(path)/$(raw.params["protocolName"])_sense_HOO_wB0_standard110_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)



