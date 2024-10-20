using HDF5, MRIReco, DelimitedFiles, BenchmarkTools, BartIO

bart = wrapper_bart(get(ENV,"TOOLBOX_PATH","/opt/software/bart-0.7.00"))
trials = parse(Int,get(ENV,"NUM_TRIALS","3"))
# change BenchmarkTools settings to match what we do in the Matlab script
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10000
BenchmarkTools.DEFAULT_PARAMETERS.samples = trials


filename = @__DIR__() * "/data/rawdata_brain_radial_96proj_12ch.h5"

rawdata = permutedims(h5read(filename, "rawdata"),[4,3,2,1]); # bart convention
traj = permutedims(h5read(filename, "trajectory"),[3,2,1]);
N = 300
toeplitz = parse(Int,get(ENV,"TOEPLITZ","1"))
oversamplingFactor = parse(Float64,get(ENV,"OVERSAMPLING","2.0"))
_,nFE,nSpokes,nCh = size(rawdata);

## Demo: NUFFT reconstruction with BART

f_sensitivity  = @__DIR__() * "/data/sensitivitiesBART.h5"

if !isfile(f_sensitivity)
  @info "Espirit"
  # inverse gridding
  img_igrid = bart(1,"nufft -i -t", traj, rawdata);
  ## Espirit coil sensitivities
  kspace_calib = bart(1,"fft -u 3", img_igrid);
  calib,emaps = bart(2,"ecalib -r 30", kspace_calib);
  sensitivity = bart(1,"slice 4 0", calib);  
  h5write(f_sensitivity, "/sensitivity", sensitivity)
else
  sensitivity = h5read(f_sensitivity, "/sensitivity")
end


## L2-SENSE reco for reference
#img_ref = @belapsed bart(1,"pics -l2 -r 0.001 -i 100 -t", traj, rawdata, smaps) evals=20


##############################
# undersampled reconstructions
##############################
@info "undersampled reco"
rf = [1,2,3,4]
img_cg = Array{ComplexF32,3}(undef,N,N,4)
times = zeros(length(rf))

for (i,d) in enumerate(rf)
  @info "r=$(d)"
  # undersample profiles
  traj_sub = traj[:,:,1:d:nSpokes]
  rawdata_sub = rawdata[:,:,1:d:nSpokes,:]

  # SENSE reconstruction while monitoring error
  img_cg[:,:,i] = bart(1,"pics -l2 -r 0.001 -i 20 -t", traj_sub, rawdata_sub, sensitivity)
  timesTrials = zeros(Float64,trials)
  for k in range(1,trials) #can't use belapsed... don't know why
    timesTrials[k] = @elapsed bart(1,"pics -l2 -r 0.001 -i 20 -t", traj_sub, rawdata_sub, sensitivity);
  end
  times[i] = minimum(timesTrials)
  @info times[i]
end

##############
# write output
##############

f_times = @__DIR__() * "/reco/recoTimes.csv"
f_img  = @__DIR__() * "/reco/images.h5"

nthreads = parse(Int,ENV["OMP_NUM_THREADS"]);
open(f_times,"a") do file
  writedlm(file, hcat("BART", nthreads, toeplitz, oversamplingFactor, transpose(times)), ',')
end

if nthreads == 1
  h5write(f_img, "/recoBART", img_cg)
end

exit()
