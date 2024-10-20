include(@__DIR__() * "/../configuration.jl")

f_times = @__DIR__() * "/reco/recoTimes.csv"
f_img  = @__DIR__() * "/reco/images.h5"
rm(f_times, force=true)
rm(f_img, force=true)
mkpath(@__DIR__() * "/reco/")

for t in threads
  @info "Run with $(t) threads"
  ENV["OMP_NUM_THREADS"] = t
  ENV["TOEPLITZ"] = 1
  ENV["OVERSAMPLING"] = 2.0
  cmd = `julia -t $t $(@__DIR__())/recoBrainBART.jl`
  @info cmd
  run(cmd)

  cmd = `julia -t $t $(@__DIR__())/recoBrainMRIReco.jl`
  @info cmd
  run(cmd)
  
  ENV["TOEPLITZ"] = 0
  ENV["OVERSAMPLING"] = 1.25 
  cmd = `julia -t $t $(@__DIR__())/recoBrainMRIReco.jl`
  @info cmd
  run(cmd)
end

include("makeImages.jl")



