using HDF5, DelimitedFiles, BenchmarkTools
using MRIReco, MRISampling, MRICoilSensitivities
using KomaMRI
trials = parse(Int,get(ENV,"NUM_TRIALS","3"))
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10000
BenchmarkTools.DEFAULT_PARAMETERS.samples = trials

filename = @__DIR__() * "/data/rawdata_brain_radial_96proj_12ch.h5"
data = permutedims(h5read(filename, "rawdata"),[3,2,1,4])
traj = permutedims(h5read(filename, "trajectory"),[3,2,1])
N = 300
Nc = 12
T = Float32
toeplitz = parse(Int,get(ENV,"TOEPLITZ","0"))
oversamplingFactor = parse(Float64,get(ENV,"OVERSAMPLING","1.25"))

@info "Threads = $(Threads.nthreads()) Toeplitz=$(toeplitz)  OversamplingFactor=$(oversamplingFactor)"   

#############################################
# load data and form Acquisition data object
##############################################
tr = Trajectory(T.(reshape(traj[1:2,:,:],2,:) ./ N), 96, 512, circular=false)
dat = Array{Array{Complex{T},2},3}(undef,1,1,1)
dat[1,1,1] = 1.e8.*reshape(data,:,12)
acqData = AcquisitionData(tr, dat, encodingSize=(N,N))

################################
# generate coil sensitivity maps
################################


f_sensitivity  = @__DIR__() * "/data/sensitivitiesMRIReco.h5"

if !isfile(f_sensitivity)
  @info "Espirit"
  acqDataCart = regrid(acqData, (N,N); cgnr_iter=3)
  sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98)
  h5write(f_sensitivity, "/sensitivity", sensitivity)
else
  sensitivity = h5read(f_sensitivity, "/sensitivity")
end

##########################
# reference reconstruction
##########################
@info "reference reco"
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 100
params[:solver] = "cgnr"
params[:toeplitz] = toeplitz == 1
params[:oversamplingFactor] = oversamplingFactor
params[:senseMaps] = Complex{T}.(reshape(sensitivity, N, N, 1, Nc))

# @time img_ref = reconstruction(acqData, params).data

##############################
# undersampled reconstructions
##############################
@info "undersampled reco"
rf = [1,2,3,4]
img_cg = Array{ComplexF32,3}(undef,N,N,4)
times = zeros(length(rf))
params[:iterations] = 20
params[:relTol] = 0.0
for i = 1:length(rf)
  @info "r=$(rf[i])"
  # undersample profiles
  global acqDataSub = convertUndersampledData(sample_kspace(acqData, T(rf[i]), "regular"))
  # SENSE reconstruction while monitoring error
  # run twice to take factor out precompilation effects
  img_cg[:,:,i] = reconstruction(acqDataSub, params).data
  #times[i] = @belapsed reconstruction(acqDataSub, params).data
  timesTrials = zeros(Float64,trials)
  for k in range(1,trials) #can't use belapsed... don't know why
    timesTrials[k] = @elapsed reconstruction(acqDataSub, params).data
  end
  times[i] = minimum(timesTrials)
  @info times[i]
end


## recon plots
width = 350
height = 300

p1 = plot_image(abs.(rotr90(img_cg[:,:,1])); title="R = 1", width=width, height=height);
p2 = plot_image(abs.(rotr90(img_cg[:,:,2])); title="R = 2", width=width, height=height);
p3 = plot_image(abs.(rotr90(img_cg[:,:,3])); title="R = 3", width=width, height=height);
p4 = plot_image(abs.(rotr90(img_cg[:,:,4])); title="R = 4", width=width, height=height);
[p1 p2 p3 p4]

## trajectory plots
tr1 = convertUndersampledData(sample_kspace(acqData, 1, "regular")).traj[1].nodes
tr2 = convertUndersampledData(sample_kspace(acqData, 2, "regular")).traj[1].nodes
tr3 = convertUndersampledData(sample_kspace(acqData, 3, "regular")).traj[1].nodes
tr4 = convertUndersampledData(sample_kspace(acqData, 4, "regular")).traj[1].nodes
s1 = scattergl(x = tr1[1,:], y=tr1[2,:],mode="markers", marker=attr(size=1, color="#EF553B"))
s2 = scattergl(x = tr2[1,:], y=tr2[2,:],mode="markers", marker=attr(size=1, color="#EF553B"))
s3 = scattergl(x = tr3[1,:], y=tr3[2,:],mode="markers", marker=attr(size=1, color="#EF553B"))
s4 = scattergl(x = tr4[1,:], y=tr4[2,:],mode="markers", marker=attr(size=1, color="#EF553B"))
p1 = plot(s1, Layout(;title="R = 1",width=400,height=400,yaxis=attr(scaleanchor="x", scaleratio=1)));
p2 = plot(s2, Layout(;title="R = 2",width=400,height=400,yaxis=attr(scaleanchor="x", scaleratio=1)));
p3 = plot(s3, Layout(;title="R = 3",width=400,height=400,yaxis=attr(scaleanchor="x", scaleratio=1)));
p4 = plot(s4, Layout(;title="R = 4",width=400,height=400,yaxis=attr(scaleanchor="x", scaleratio=1)));
[p1 p2 p3 p4]

## sensitivity plots
s = [sensitivity[:,:,1,i] for i=1:12]
s = [s[1] s[2] s[3] s[4]; s[5] s[6] s[7] s[8]; s[9] s[10] s[11] s[12]]
plot_image(abs.(rotr90(s)); title="sensitivity - espirit")




params = Dict{Symbol,Any}()
params[:reconSize] = (N, N)
params[:densityWeighting] = true

image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]

rf = [1,2,3,4]
img_cg = Array{ComplexF32,3}(undef,N,N,4)
for i = 1:length(rf)
  @info "r=$(rf[i])"
  # undersample profiles
  global acqDataSub = convertUndersampledData(sample_kspace(acqData, T(rf[i]), "regular"))
  # SENSE reconstruction while monitoring error
  # run twice to take factor out precompilation effects
  img_cg[:,:,i] = reconstruction(acqDataSub, params).data[:,:,1,1,1,1]
end
p1 = plot_image(abs.(rotr90(img_cg[:,:,1])); title="R = 1", width=width, height=height);
p2 = plot_image(abs.(rotr90(img_cg[:,:,2])); title="R = 2", width=width, height=height);
p3 = plot_image(abs.(rotr90(img_cg[:,:,3])); title="R = 3", width=width, height=height);
p4 = plot_image(abs.(rotr90(img_cg[:,:,4])); title="R = 4", width=width, height=height);
[p1 p2 p3 p4]
##############
# write output
##############

f_times = @__DIR__() * "/reco/recoTimes.csv"
f_img  = @__DIR__() * "/reco/images.h5"

open(f_times,"a") do file
  writedlm(file, hcat("MRIReco", Threads.nthreads(), toeplitz, oversamplingFactor, transpose(times)), ',')
end

if Threads.nthreads() == 1
  h5write(f_img, "/recoMRIReco$(toeplitz)", img_cg)
end

exit()
