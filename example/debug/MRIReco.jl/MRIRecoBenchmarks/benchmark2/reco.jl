using HDF5, BartIO, DelimitedFiles, BenchmarkTools, ProfileView
using ImageQualityIndexes:assess_ssim
using MRIReco, MRISampling, MRICoilSensitivities, ImageUtils

bart = wrapper_bart(get(ENV,"TOOLBOX_PATH","/opt/software/bart-0.7.00"))
trials = parse(Int,get(ENV,"NUM_TRIALS","1"))
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10000
BenchmarkTools.DEFAULT_PARAMETERS.samples = trials

T = Float32

phantom = bart(1,"phantom");
phantom = bart(1,"noise -n0.005",phantom)
N = 128
NCh = 8
phant3D = bart(1,"phantom -3 -x$N -s$NCh");
phant3D = bart(1,"noise -n0.005",phant3D)
phant3D_rss = bart(1,"rss 8",phant3D)
kbart = bart(1,"fft -u 7",phant3D);
mask = bart(1,"poisson -Y $N -Z $N -y1.2 -z1.2 -C 20 -v -V5");
kbart_u = kbart .* mask;

################################
# generate coil sensitivity maps
################################


f_sensitivity  = @__DIR__() * "/data/sensitivities.h5"

if !isfile(f_sensitivity)
  @info "ecalib"
  sensitivity =  bart(1,"ecalib -m1 -c0", kbart_u);
  h5write(f_sensitivity, "/sensitivity", sensitivity)
else
  sensitivity = h5read(f_sensitivity, "/sensitivity")
end

imFully = bart(1,"pics -d5 -i1 l2 -r0",kbart,sensitivity);

################################
# reconstruction BART
################################

imBART = bart(1,"pics -d5 -i30 -RW:7:0:0.01", kbart_u, sensitivity)
timesTrials = zeros(Float64,trials)
for k in range(1,trials) 
    timesTrials[k] = @elapsed bart(1,"pics -d5 -i30 -RW:7:0:0.01", kbart_u, sensitivity)
end
timeBART = minimum(timesTrials)
@info timeBART

RMSE_bart = MRIReco.norm(vec(abs.(imBART))-vec(abs.(imFully)))/MRIReco.norm(vec(abs.(imFully)))
ssim_bart = round(assess_ssim(abs.(imBART[:,:,80]),abs.(imFully[:,:,80])),digits=3)


################################
# reconstruction MRIReco.jl
################################

tr = MRIBase.CartesianTrajectory3D(T, N, N, numSlices=N, TE=T(0), AQ=T(0))
kdata_j = [reshape(Complex{T}.(kbart),:,NCh) for i=1:1, j=1:1, k=1:1]
acq = AcquisitionData(tr, kdata_j, encodingSize=(N,N,N))

# Do we need this ????
#params = Dict{Symbol, Any}()
#params[:reco] = "direct"
#params[:reconSize] = tuple(acq.encodingSize...)
#Ireco = reconstruction(acq, params)
#Isos = mergeChannels(Ireco)
#heatmap(abs.(Isos[:,:,80]), c=:grays, aspect_ratio = 1,legend = :none , axis=nothing,showaxis = false)

# find indices
I = findall(x->x==1,abs.(repeat(mask,N,1,1)))
subsampleInd = LinearIndices((N,N,N))[I]

acqCS = deepcopy(acq);
acqCS.subsampleIndices[1] = subsampleInd
acqCS.kdata[1,1,1] = acqCS.kdata[1,1,1][subsampleInd,:]


params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:senseMaps] = Complex{T}.(sensitivity);

params[:solver] = "fista"
params[:sparseTrafoName] = "Wavelet"
params[:regularization] = "L1"
params[:λ] = T(0.01) # 5.e-2
params[:iterations] = 30
params[:normalize_ρ] = false
params[:ρ] = T(0.95)
#params[:relTol] = 0.1
params[:normalizeReg] = true


imMRIReco = reconstruction(acqCS, params).data;

for k in range(1,trials) 
    timesTrials[k] = @elapsed reconstruction(acqCS, params);
end
timeMRIReco = minimum(timesTrials)
@info timeMRIReco



##############
# write output
##############

f_times = @__DIR__() * "/reco/recoTimes.csv"
f_img  = @__DIR__() * "/reco/images.h5"
nthreads = parse(Int,ENV["OMP_NUM_THREADS"]);

open(f_times,"a") do file
  writedlm(file, hcat("BART", nthreads, timeBART), ',')
  writedlm(file, hcat("MRIReco", Threads.nthreads(), timeMRIReco), ',')
end

if Threads.nthreads() == 1
  h5write(f_img, "/imFully", imFully)
  h5write(f_img, "/recoBART", imBART)
  h5write(f_img, "/recoMRIReco", imMRIReco)
end

exit()
