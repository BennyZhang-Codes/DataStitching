using KomaHighOrder, MRIReco, MRICoilSensitivities, PyPlot, MAT
import KomaHighOrder.MRIBase: AcquisitionData

#=
prepare data for synchronization
1. Multi-Echo GRE (ME-GRE) data: 
    (1) Coil-Sensitivity Map (CSM)
    (2) ΔB₀ map
2. Spiral imaging datas 
    (1) MRI Signal data

Notes: directory structure:
    - path/seq/*.seq     # sequence files
    - path/mrd/*.mrd     # MRD files
    - path/syn/*.mat     # output directory
=#

path         = "/home/jyzhang/Desktop/pulseq/20241010_skope_fa90/invivo"
gre_seq_file = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("gres6", f)][1]
gre_mrd_file = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin("gres6", f)][1]

if ispath("$(path)/syn") == false mkpath("$(path)/syn") end

T              = Float32
kdata, ktraj, kdims, shape, TE, ReadoutMode, acqData, raw = read_gre(gre_seq_file, gre_mrd_file);
nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape;

# IFFT of gre
imgs = convert_ifft(kdata, dims=[2,3]);
imgs = CoilCombineSOS(abs.(imgs), 1);
imgs = permutedims(imgs, [3,1,2]);
fig = plt_images(imgs,width=7.5, height=5, vminp=0, vmaxp=99)

############
# espirit
############
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.99);  # (nX, nY, 1, nCha)
csm = permutedims(sensitivity, [2,1,4,3])[:,:,:,1];# (nX, nY, 1, nCha) => (nY, nX, nCha, 1) => (nY, nX, nCha)
fig = plt_images(permutedims(abs.(csm), [3,1,2]),width=10, height=5)  # (nCha, nY, nX)

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
fig = plt_B0map(b0, width=5, height=4)


MAT.matwrite("$(path)/syn/syn_$(raw.params["protocolName"]).mat", Dict("csm" => csm, "b0" => b0, "mask" => mask))


extract_data(r::Int64, path::String) = begin
    seq_file     = "$(path)/seq/" * [f for f in readdir("$(path)/seq") if occursin("r$(r)", f)][1]
    mrd_file     = "$(path)/mrd/" * [f for f in readdir("$(path)/mrd") if occursin("r$(r)", f)][1]

    seq = read_seq(seq_file)[end-9:end-3]; 
    # times = KomaMRIBase.get_adc_sampling_times(seq);
    FOV        = seq.DEF["FOV"]
    matrixSize = Int.(seq.DEF["matrixSize"])

    raw = RawAcquisitionData(ISMRMRDFile(mrd_file))
    raw.params["trajectory"] = "custom";
    shape = get_ksize(raw);
    nCha, nZ, nY, nX, nAvg, nSli, nCon, nPha, nRep, nSet, nSeg = shape; println(shape);
    kdata = Complex{T}.(get_kdata(raw, shape));
    kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...));
    kdims = [mrddims[idx] for idx in 1:length(shape) if shape[idx]>1];

    data = permutedims(reshape(kdata, nCha, nX*nSet), [2,1]);

    MAT.matwrite("$(path)/syn/syn_$(raw.params["protocolName"]).mat", Dict("data" => data, "FOV" => FOV, "matrixSize" => matrixSize))
end
for r = 1:4
    extract_data(r, path)
end