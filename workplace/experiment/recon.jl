using KomaHighOrder, MRIReco, PlotlyJS
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory, subsampleIndices, rawdata
import KomaHighOrder.MRIBase: uniqueidx, encSteps1, encSteps2


path       = "/home/jyzhang/Desktop/pulseq/20240624/"
seq_file   = path * "gre_fov150_150.seq"
raw_file   = path * "meas_MID00309_FID40507_pulseqt_gre_fov150_150_coilsens.mrd"
name       = "gre_fov150_150_coilsens"
T          = Float32
N          = 150
# get signal data
raw              = RawAcquisitionData(ISMRMRDFile(raw_file))
# raw = demo_raw(BlochHighOrder("000"))   ##
NumberOfSamples  = Int64(raw.profiles[1].head.number_of_samples);
NumberOfProfiles = Int64(length(raw.profiles));
NumberOfChannels = raw.params["receiverChannels"];
raw.params["trajectory"] = "custom";
kdata = rawdata(raw);

# get trajectory data
seq          = read_seq(seq_file);
# seq = demo_seq()                      ##

_, ktraj_adc = get_kspace(seq);
tr           = Trajectory(T.(ktraj_adc[:,1:2]'), NumberOfProfiles, NumberOfSamples, circular=false);
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
savefig(p_sos, path * "$(name)_sos.svg", width=450, height=420, format="svg")

# all  channels
p_cha = plot_imgs_subplots(images, 4, 8; title="All channels, $(name)", height=400, width=800)
savefig(p_cha, path * "$(name)_32cha.svg", width=800, height=400, format="svg")

# single channel
# plot_image(images[:,:,1]; )

# s = scatter(x=ktraj_adc'[1,:], y=ktraj_adc'[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)