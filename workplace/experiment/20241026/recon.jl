using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

T = Float64

path         = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo"
r=4
seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("r$(r)", f)][1]
mrd_file     = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin("r$(r)", f)][1]
syn_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin("r$(r)", f)][1]
gre_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(r"^syn.*gre.*mat$", f)][1]

outpath = "$(path)/out"; if ispath(outpath) == false mkpath(outpath) end

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
kStandard  = matread(syn_file)["ksphaStandard_syn"]/2π;
kStitched  = matread(syn_file)["ksphaStitched_syn"]/2π;

nY, nX, nCha = size(csm);
# csm = permutedims(csm, [2,1,3]); # (nX, nY, nCha)
# csm = Complex{T}.(reshape(csm[end:-1:1,:,:,:], nX, nY, 1, nCha)); # reverse the x-axis


nSample, nCha = size(data);
tr           = Trajectory(T.(k_adc[:,1:2]'), 1, nSample, circular=false, times=times);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = data;
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));
C = maximum(2*abs.(acqData.traj[1].nodes[:]));  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C;


tr_nominal  = Trajectory(k_adc'[1:3,:],   1, nSample; circular=false, times=times);
tr_Stitched = Trajectory(kStitched'[:,:], 1, nSample; circular=false, times=times);
tr_Standard = Trajectory(kStandard'[:,:], 1, nSample; circular=false, times=times);
########################################################################
# SENSE
########################################################################
# solver = "cgnr"
# reg = "L2"
# iter = 20
# λ = 1e-3
# recParams = Dict{Symbol,Any}()
# recParams[:reconSize]        = (nX, nY)
# recParams[:densityWeighting] = true
# recParams[:reco] = "multiCoil"
# recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
# recParams[:λ] = λ
# recParams[:iterations] = iter
# recParams[:solver] = solver
# recParams[:oversamplingFactor] = 2
# csmap = permutedims(csm, [2,1,3]); # (nX, nY, nCha)
# csmap = Complex{T}.(reshape(csmap[end:-1:1,:,:,:], nX, nY, 1, nCha)); # reverse the x-axis
# recParams[:senseMaps] = csmap;

# @time rec = reconstruction(acqData, recParams).data[:,:];
# fig = plt_image(rotr90(abs.(rec)))

# fig.savefig("$(path)/$(raw.params["protocolName"])_sense_$(solver)_$(reg)_$(iter)_$(λ).png", dpi=300, bbox_inches="tight", pad_inches=0)



########################################################################
# HighOrderOp
########################################################################
b0map = Matrix(b0map); 
Nblocks=6;
grid=5;

solver = "cgnr"
reg = "L2"
iter = 2
λ = 1e-2
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (nX, nY)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
# recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = Complex{T}.(reshape(csm, nX, nY, 1, nCha));
S = SensitivityOp(reshape(Complex{T}.(csm),:,nCha),1);

# tr_Stitched = Trajectory(kStitched'[:,:].*(1), 1, nSample; circular=false, times=times);
Op0 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched, BlochHighOrder("010"); fieldmap=b0map, Nblocks=Nblocks, grid=grid, verbose=true);
recParams[:encodingOps] = reshape([DiagOp(Op0, nCha) ∘ S], 1,1);
@time rec = reconstruction(acqData, recParams).data[:,:];
plt_image(abs.(rec))



# 1. nominal trajectory, BlochHighOrder("000")
Op1 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);
# 2. stitched trajectory, BlochHighOrder("011")
Op2 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("011"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);
# 3. stitched trajectory, BlochHighOrder("110")
Op3 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);
# 4. stitched trajectory, BlochHighOrder("111")
Op4 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);
# 5. standard trajectory, BlochHighOrder("011")
Op5 = HighOrderOp((nX, nY), tr_nominal, tr_Standard , BlochHighOrder("011"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);
# 6. standard trajectory, BlochHighOrder("110")
Op6 = HighOrderOp((nX, nY), tr_nominal, tr_Standard , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);
# 7. standard trajectory, BlochHighOrder("111")
Op7 = HighOrderOp((nX, nY), tr_nominal, tr_Standard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=b0map, grid=grid, verbose=false);

# 8. nominal trajectory, BlochHighOrder("000")
Op8 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("000"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
# 9. stitched trajectory, BlochHighOrder("011")
Op9 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("011"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
# 10. stitched trajectory, BlochHighOrder("110")
Op10 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
# 11. stitched trajectory, BlochHighOrder("111")
Op11 = HighOrderOp((nX, nY), tr_nominal, tr_Stitched , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
# 12. standard trajectory, BlochHighOrder("011")
Op12 = HighOrderOp((nX, nY), tr_nominal, tr_Standard , BlochHighOrder("011"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
# 13. standard trajectory, BlochHighOrder("110")
Op13 = HighOrderOp((nX, nY), tr_nominal, tr_Standard , BlochHighOrder("110"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
# 14. standard trajectory, BlochHighOrder("111")
Op14 = HighOrderOp((nX, nY), tr_nominal, tr_Standard , BlochHighOrder("111"); Nblocks=Nblocks, fieldmap=b0map.*0, grid=grid, verbose=false);
Ops = [Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8, Op9, Op10, Op11, Op12, Op13, Op14];

imgs = Array{Complex{T},3}(undef, length(Ops), nX, nY);
labels = [
            "wB0_nominal",
            "wB0_stitched_011",
            "wB0_stitched_110",
            "wB0_stitched_111",
            "wB0_standard_011",
            "wB0_standard_110",
            "wB0_standard_111",
            "woB0_nominal",
            "woB0_stitched_011",
            "woB0_stitched_110",
            "woB0_stitched_111",
            "woB0_standard_011",
            "woB0_standard_110",
            "woB0_standard_111",
          ];
for idx in eachindex(Ops)
    recParams[:encodingOps] = reshape([DiagOp(Ops[idx], nCha) ∘ S], 1,1);
    @time rec = reconstruction(acqData, recParams).data[:,:];
    imgs[idx, :, :] = rec;
    fig = plt_image(abs.(rec); title=labels[idx])
end

MAT.matwrite("$(path)/out/r$(r)_$(solver)_$(reg)_$(λ)_$(iter).mat", Dict("imgs"=>imgs, "labels"=>labels))








