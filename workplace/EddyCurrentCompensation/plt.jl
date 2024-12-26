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

labels = ["model", "Stitched", "Standard"]
xs = [times.*1e3]
ys = [vec(k0_ecc), kStitched[:, 1]*2π, kStandard[:, 1]*2π]
fig = plt_plot(ys; labels=labels, xs=xs, ylabel=L"k_0 \ [rad]", xlabel=L"Time \ [ms]", width=16, height=8, fontsize_label=10, fontsize_legend=10, fontsize_ticklabel=8)
fig.savefig("/home/jyzhang/Desktop/ECC_k0.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)
