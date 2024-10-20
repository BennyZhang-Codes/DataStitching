using KomaHighOrder, MRIReco, PlotlyJS
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory, subsampleIndices, rawdata
import KomaHighOrder.MRIBase: uniqueidx, encSteps1, encSteps2
using MRIReco, MRICoilSensitivities

path       = "$(@__DIR__)/edit_phantom/coilsensitivity/"
seq_file   = path * "gre_fov150_150.seq"
raw_file   = path * "meas_MID00309_FID40507_pulseqt_gre_fov150_150_coilsens.mrd"
name       = "gre_fov150_150_coilsens"
T          = Float32
N          = 150
# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(raw_file))

NumberOfSamples  = Int64(raw.profiles[1].head.number_of_samples);
NumberOfProfiles = Int64(length(raw.profiles));
NumberOfChannels = raw.params["receiverChannels"];
raw.params["trajectory"] = "custom";
kdata = rawdata(raw);

# get trajectory data
seq          = read_seq(seq_file);

_, ktraj_adc = get_kspace(seq);
tr           = Trajectory(T.(ktraj_adc[:,1:2]'), NumberOfProfiles, NumberOfSamples, circular=false, cartesian=true);
dat          = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1]   = reshape(kdata,:,NumberOfChannels);
acqData      = AcquisitionData(tr, dat, encodingSize=(N,N));

# acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
Nx = Ny = N
recParams = Dict{Symbol,Any}()
recParams[:reconSize]        = (Nx, Ny)
recParams[:densityWeighting] = true
rec      = reconstruction(acqData, recParams);
images   = reshape(rec.data, Nx, Ny, :);
images   = (abs.(images) * prod(size(images)[1:2]));

# Sum-Of-Square
p_sos = plot_img(sqrt.(sum(images.^2; dims=3))[:,:,1]; title="SOS, $(name)", width=450, height=420)
# savefig(p_sos, path * "$(raw.params["protocolName"])_sos.svg", width=450, height=420, format="svg")

# all  channels
p_cha = plot_imgs_subplots(images, 4, 8; title="All channels, $(name)", height=400, width=800)
# savefig(p_cha, path * "$(raw.params["protocolName"])_32cha.svg", width=800, height=400, format="svg")

# single channel
# plot_image(images[:,:,1]; )

# s = scatter(x=ktraj_adc'[1,:], y=ktraj_adc'[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)


# espirit
# acqDataCart = regrid(acqData, (Nx, Ny); cgnr_iter=3);
sensitivity = espirit(acqData, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98);
plot_imgs_subplots(real(sensitivity[:,:,1,:]), 4, 8)
plot_img(abs.(sqrt.(sum(sensitivity[:,:,1,:].^2; dims=3))[:,:,1]))
p_smap_espirit = plot_imgs_subplots(abs.(sensitivity[:,:,1,:]), 4, 8; title="Coil Sensitivity (espirit)")
savefig(p_smap_espirit,  path * "$(raw.params["protocolName"])-CoilSens_espirit.svg", format="svg", height=400, width=800)

s = sensitivity[:,:,1,:]
MAT.matwrite(path * "coilsensmap_32cha.mat",  Dict("coilsensmap_32cha"=>s))


import ImageTransformations: imresize


n = 250
s_real = real(s)
s_imag = imag(s)
ss = imresize(s_real, (n,n,32)) .+ imresize(s_imag, (n,n,32)) * im

plot_imgs_subplots(abs.(s1), 4, 8)


function RealCoilSensitivity_32cha(Nx::Int64, Ny::Int64)
    sensitivity = MAT.matread("$(@__DIR__)/coilsensmap_32cha.mat")["coilsensmap_32cha"]
    out = imresize(sensitivity, (Nx, Ny, 32))
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
    out = out./ norm
    return out
end

plot_img(abs.(sqrt.(sum(smap.^2; dims=3))[:,:,1]))
norm = sqrt.(sum(abs.(smap) .^ 2, dims=3))
smap_norm = smap./ norm
plot_img(abs.(sqrt.(sum(smap_norm.^2; dims=3))[:,:,1]))

n=1000
smap = RealCoilSensitivity_32cha(n, n)
plot_imgs_subplots(abs.(smap), 4, 8)