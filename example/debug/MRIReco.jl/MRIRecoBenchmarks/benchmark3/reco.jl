using HDF5, BartIO, DelimitedFiles, BenchmarkTools, ProfileView
using ImageQualityIndexes:assess_ssim
using MRIReco, MRISampling, MRICoilSensitivities, ImageUtils

bart = wrapper_bart(get(ENV,"TOOLBOX_PATH","/opt/software/bart-0.7.00"))
trials = parse(Int,get(ENV,"NUM_TRIALS","1"))
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10000
BenchmarkTools.DEFAULT_PARAMETERS.samples = trials

T = Float32
N = 128
NCh = 8

sensitivity = bart(1,"phantom -3 -x$N -S$NCh")
phant3D = bart(1,"phantom -3 -x$N -s$NCh")
phant3D = bart(1,"noise -n0.005", phant3D)
phant3D_rss = bart(1,"rss 8", phant3D)
kbart = bart(1,"fft -u 7", phant3D)


################################
# espirit BART
################################

sensitivityBART = bart(1,"ecalib -m1 -c0 -k 6 -r24", kbart)
timesTrials = zeros(Float64,trials)
for k in range(1,trials) 
    timesTrials[k] = @elapsed bart(1,"ecalib -m1 -c0 -k 6 -r24", kbart)
end
timeBART = minimum(timesTrials)
@info timeBART

################################
# espirit MRIReco.jl
################################

"""
Documentation : Crop the central area for 4D array
"""
function crop(A::Array{T,4}, s::NTuple{3,Int64}) where {T}
    nx, ny, nz = size(A)
    idx_x = div(nx, 2)-div(s[1], 2)+1:div(nx, 2)-div(s[1], 2)+s[1]
    idx_y = div(ny, 2)-div(s[2], 2)+1:div(ny, 2)-div(s[2], 2)+s[2]
    idx_z = div(nz, 2)-div(s[3], 2)+1:div(nz, 2)-div(s[3], 2)+s[3]
    return A[idx_x, idx_y, idx_z,:]
end

sensitivityMRIReco = espirit(crop(kbart, (24,24,24)), (N,N,N), (6,6,6), nmaps=1, eigThresh_2=0.0)

for k in range(1,trials) 
    timesTrials[k] = @elapsed espirit(crop(kbart, (24,24,24)), (N,N,N), (6,6,6), nmaps=1, eigThresh_2=0.0)
end
timeMRIReco = minimum(timesTrials)
@info timeMRIReco



##############
# write output
##############

f_times = @__DIR__() * "/reco/recoTimes.csv"
f_img  = @__DIR__() * "/reco/images.h5"
nthreads = parse(Int,get(ENV,"OMP_NUM_THREADS","1")) 

open(f_times,"a") do file
  writedlm(file, hcat("BART", nthreads, timeBART), ',')
  writedlm(file, hcat("MRIReco", Threads.nthreads(), timeMRIReco), ',')
end

if Threads.nthreads() == 1
  h5write(f_img, "/sensitivity", sensitivity)
  h5write(f_img, "/sensitivityBART", sensitivityBART)
  h5write(f_img, "/sensitivityMRIReco", sensitivityMRIReco)
end

exit()
