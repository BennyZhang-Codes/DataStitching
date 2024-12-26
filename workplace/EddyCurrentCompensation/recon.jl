using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

using CUDA

idx = 1;
CUDA.device!(1)


T = Float64
path         = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/"

seqs = ["7T_1mm-200-r4_max51-fa90.seq",
        "7T_1mm-200-r3_max51-fa90.seq",
        "7T_1mm-200-r2_max51-fa90.seq",
        "7T_0.71mm-280-r2_max51-fa90.seq",
        "7T_0.6mm-332-r3_max51-fa90.seq",
        "7T_0.5mm-400-r4_max51-fa90.seq"]
mrds = ["r4_1p0_standard", "r3_1p0_standard", "r2_1p0_standard", "r2_0p71_standard", "r3_0p6_standard", "r4_0p5_standard"]
# gres = ["gres6_1p0_standard", "gres6_1p0_standard", "gres6_1p0_standard", "gres6_1p0_standard", "gres6_1p0_standard", "gres6_1p0_standard"];

mris = ["r4_1p0_standard", "r3_1p0_standard", "r2_1p0_standard", "r2_0p71_standard", "r3_0p6_standard", "r4_0p5_standard"];

seq_pro = mris[idx]

seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin(seqs[idx], f)][1]
mrd_file     = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin(mrds[idx], f)][1]
syn_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(mris[idx], f)][1]
gre_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin("gres6_1p0_standard.mat", f)][1] 
ECC_file     = "$(path)/ecc/" * [f for f in readdir("$(path)/ecc") if occursin(mris[idx]*".mat", f)][1] 
# syn_rep1_file     = "$(path)/syn/" * [f for f in readdir("$(path)/syn") if occursin(mris[idx]*"_rep1", f)][1]

outpath = "$(@__DIR__)/workplace/EddyCurrentCompensation/out"; if ispath(outpath) == false mkpath(outpath) end

seq = read_seq(seq_file)[end-9:end-3]; 
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
TE = seq.DEF["TE"];
_, k_adc = get_kspace(seq);
times = KomaMRIBase.get_adc_sampling_times(seq);
times = times .- times[1] .+ TE;

# Coil-Sensitivity Map (CSM), ΔB0 map, mask
b0         = matread(gre_file)["b0"];
csm        = matread(gre_file)["csm"];    # (nY, nX, nCha)
mask       = matread(gre_file)["mask"];

# MRI signal data, FOV, matrix size, synchronized dynamic fields (kspha, rad, rad/m⁻¹, rad/m⁻²)
data       = matread(syn_file)["data"];
FOV        = matread(syn_file)["FOV"];
matrixSize = matread(syn_file)["matrixSize"];
kStandard  = -matread(syn_file)["ksphaStandard_syn"]/2π;
kStitched  = -matread(syn_file)["ksphaStitched_syn"]/2π;

nY, nX, nZ = matrixSize;
nSample, nCha = size(data);

b0 = imresize_real(b0, (nY, nX));
# csm = imresize_complex(csm, (nY, nX, nCha));
# norm = sqrt.(sum(abs.(csm) .^ 2, dims=3));
# csm = csm./ norm;
# csm[isnan.(csm)] .= 0 + 0im;

csm = imresize(csm, (nY, nX, nCha))


# csm = permutedims(csm, [2,1,3]); # (nX, nY, nCha)
# csm = Complex{T}.(csm[end:-1:1,:,:]); # reverse the x-axis

k0_ecc = matread(ECC_file)["phase_drift"];

fig, ax = plt.subplots(1,1)
ax.plot(k0_ecc', label="model")
ax.plot(kStitched[:, 1]*2π, label="Stitched")
ax.legend()

nSample, nCha = size(data);
tr           = Trajectory(T.(k_adc'[1:2, :]), 1, nSample, circular=false, times=times);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);

# data = data .* exp.(-1im.*k0_ecc)';
# data = data .* exp.(-1im.*kStitched[:, 1]*2π);

dat[1,1,1]   = data;
acqData      = AcquisitionData(tr, dat, encodingSize=(nX, nY));
C = maximum(2*abs.(acqData.traj[1].nodes[:]));  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C;


tr_nominal  = Trajectory(-k_adc'[1:3,:],  1, nSample; circular=false, times=times);
tr_Stitched = Trajectory(kStitched'[:,:], 1, nSample; circular=false, times=times);
tr_Standard = Trajectory(kStandard'[:,:], 1, nSample; circular=false, times=times);

########################################################################
# HighOrderOp
########################################################################
Δx = Δy = FOV[1]/matrixSize[1]


b0 = Matrix(b0); 
Nblocks=40;
grid=3;  # 5
verbose = true;

solver = "cgnr"; reg = "L2"; iter = 10; λ = 1e-3
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (nX, nY)
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:senseMaps] = Complex{T}.(reshape(csm, nX, nY, 1, nCha));
# S = SensitivityOp(reshape(Complex{T}.(csm),:,nCha),1);

BHO = BlochHighOrder("000")
Op_i2 = HighOrderOp_i2((nX, nY), tr_nominal, tr_Stitched, BHO; Δx=Δx, Δy=Δy, 
                        Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=b0.*0, grid=grid, use_gpu=true, verbose=verbose);
@time x = recon_HOOp(Op_i2, acqData, recParams);
plt_image(abs.(x), title=BHO.name, color_facecolor="#000000", color_label="#FFFFFF", vmaxp=99.9, vmin=0)
plt_image(angle.(x), title=BHO.name, color_facecolor="#000000", color_label="#FFFFFF", vmax=π, vmin=-π)

# fig = plt.subplots(1,1)
# plot(k_adc[:, 1], k_adc[:,2])
# plot(kStitched[:, 2], kStitched[:,3])
# plot(kStandard[:, 2], kStandard[:,3])

labels = ["000", "100", "010", "001", "110", "101", "011", "111"]
imgs = Array{Complex{T},3}(undef, length(labels), nX, nY);
for idx_l in eachindex(labels)
    BHO = labels[idx_l]
    Op = HighOrderOp_i2((nX, nY), tr_nominal, tr_Stitched, BlochHighOrder(BHO); Δx=Δx, Δy=Δy, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=b0, grid=grid, use_gpu=true, verbose=verbose);
    @time x = recon_HOOp(Op, acqData, recParams);
    imgs[idx_l, :, :] = x;
    fig_mag = plt_image(  abs.(x); title=BHO, color_facecolor="#000000", color_label="#FFFFFF", vmaxp=99.9, vmin=0)
    fig_mag.savefig("$(outpath)/$(seq_pro)_$(BHO)_mag.png", dpi=600, transparent=false, bbox_inches="tight", pad_inches=0)
    fig_pha = plt_image(angle.(x); title=BHO, color_facecolor="#000000", color_label="#FFFFFF", vmax=π, vmin=-π)
    fig_pha.savefig("$(outpath)/$(seq_pro)_$(BHO)_pha.png", dpi=600, transparent=false, bbox_inches="tight", pad_inches=0)
end
MAT.matwrite("$(outpath)/$(seq_pro)_$(solver)_$(reg)_$(λ)_$(iter).mat", Dict("imgs"=>imgs, "labels"=>labels))

