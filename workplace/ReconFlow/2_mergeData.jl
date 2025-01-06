import KomaHighOrder: get_grads
using MAT
T = Float32

path         = "/home/jyzhang/Desktop/pulseq/20241104_ABDL/"
out_path     = "$(path)/data"
if ispath(out_path) == false mkpath(out_path) end

seqs = ["7T_1mm-200-r4_max51-fa90.seq",
        "7T_1mm-200-r3_max51-fa90.seq",
        "7T_1mm-200-r2_max51-fa90.seq",
        "7T_0.71mm-280-r2_max51-fa90.seq",
        "7T_0.6mm-332-r3_max51-fa90.seq",
        "7T_0.5mm-400-r4_max51-fa90.seq"]

mrds = ["meas_MID00115_FID53003_pulseq_v0_r4_1p0_standard.mrd",
        "meas_MID00114_FID53002_pulseq_v0_r3_1p0_standard.mrd",
        "meas_MID00113_FID53001_pulseq_v0_r2_1p0_standard.mrd",
        "meas_MID00118_FID53006_pulseq_v0_r2_0p71_standard.mrd",
        "meas_MID00119_FID53007_pulseq_v0_r3_0p6_standard.mrd",
        "meas_MID00126_FID53014_pulseq_v0_r4_0p5_standard.mrd"]

gres = ["gre_1p0_200_r1.mat",
        "gre_1p0_200_r1.mat",
        "gre_1p0_200_r1.mat",
        "gre_0p71_280_r1.mat",
        "gre_0p6_332_r1.mat",
        "gre_0p5_400_r1.mat"]

dfcs = ["7T_1p0_200_r4.mat",
        "7T_1p0_200_r3.mat",
        "7T_1p0_200_r2.mat",
        "7T_0p71_280_r2.mat",
        "7T_0p6_332_r3.mat",
        "7T_0p5_400_r4.mat"]

filenames = ["1p0_200_r4",
             "1p0_200_r3",
             "1p0_200_r2",
             "0p71_280_r2",
             "0p6_332_r3",
             "0p5_400_r4"]


for idx in eachindex(seqs)
# idx = 1
seq = seqs[idx]
mrd = mrds[idx]
gre = gres[idx]
dfc = dfcs[idx]
filename = filenames[idx]

seq_file = "$(path)/seq/$(seq)" 
mrd_file = "$(path)/mrd/$(mrd)"
gre_file = "$(path)/gre/$(gre)"
dfc_file = "$(path)/dfc/$(dfc)"

@info "filename: $(filename)"
@info "seq file: $(seq_file)"
@info "mrd file: $(mrd_file)"
@info "gre file: $(gre_file)"
@info "dfc file: $(dfc_file)"

seq = read_seq(seq_file)[end-9:end-3];
seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
seq_dt      = seq.DEF["adc_DwellTime"];
FOV         = seq.DEF["FOV"]
matrixSize  = Int.(seq.DEF["matrixSize"])
TE          = seq.DEF["TE"];

times = KomaMRIBase.get_adc_sampling_times(seq);
times = times .- times[1] .+ TE;

raw = RawAcquisitionData(ISMRMRDFile(mrd_file))
raw.params["trajectory"] = "custom";
shape = get_ksize(raw);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape);
kdata = Complex{T}.(get_kdata(raw, shape));

# kdata = mean(kdata, dims=9);
kdata = kdata[:,:,:,:,:,:,:,:,1,:,:];

kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];

kdata = permutedims(reshape(kdata, nCha, nX*nSet), [2,1]);

nSample = nX*nSet;
datatime = seq_dt .* collect(0:nSample-1);

# gre data
gre_data = matread(gre_file);
gre_img  = gre_data["img"];
gre_csm  = gre_data["csm"];
gre_b0   = gre_data["b0"];
gre_mask = gre_data["mask"];

# dynamic field data
dfc_data          = matread(dfc_file);
dfc_dt            = dfc_data["dt"];
dfc_ksphaStitched = dfc_data["ksphaStitched"]./2π;  # rad to Hz
dfc_ksphaStandard = dfc_data["ksphaStandard"]./2π;
dfc_startStitched = dfc_data["delayStitched"];
dfc_startStandard = dfc_data["delayStandard"];
dfc_nSampleAllSegStitched = dfc_data["nSampleAllSegStitched"];
dfc_nSampleAllSegStandard = dfc_data["nSampleAllSegStandard"];

# get nominal trajectory
k = nothing
for s in seq
    if KomaHighOrder.is_ADC_on(s)
        t = collect(0:dfc_dt:s.DUR[1]);
        gx, gy, gz = get_grads(s, t)
        kx = cumtrapz(ones(length(t)-1)'*dfc_dt, gx')
        ky = cumtrapz(ones(length(t)-1)'*dfc_dt, gy')
        kz = cumtrapz(ones(length(t)-1)'*dfc_dt, gz')
        k = [kx' ky' kz'] * γ
    end
end

ksphaNominal = zeros(Float64, size(k, 1), 9);
ksphaNominal[:,2:4] = k[:,1:3];
startNominal = 0.;

MAT.matwrite("$(out_path)/$(filename).mat", 
    Dict("kdata" => kdata, "datatime" => datatime, "FOV" => FOV, "matrixSize" => matrixSize, "seq_dt" => seq_dt, "TE" => TE,
        "ksphaNominal" => ksphaNominal, "startNominal" => startNominal,
        "gre_img" => gre_img, "gre_csm" => gre_csm, "gre_b0" => gre_b0, "gre_mask" => gre_mask,
        "dfc_dt" => dfc_dt, "dfc_ksphaStitched" => dfc_ksphaStitched, "dfc_ksphaStandard" => dfc_ksphaStandard,
        "dfc_startStitched" => dfc_startStitched, "dfc_startStandard" => dfc_startStandard,
        "dfc_nSampleAllSegStitched" => dfc_nSampleAllSegStitched, "dfc_nSampleAllSegStandard" => dfc_nSampleAllSegStandard))

end



