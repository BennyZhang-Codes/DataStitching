include(@__DIR__() * "/../configuration.jl")


f_times = @__DIR__() * "/reco/recoTimes.csv"
f_img  = @__DIR__() * "/reco/images.h5"
rm(f_times, force=true)
rm(f_img, force=true)
mkpath(@__DIR__() * "/reco/")
mkpath(@__DIR__() * "/data/")

for t in threads
  @info "Run with $(t) threads"
  ENV["OMP_NUM_THREADS"] = t
  cmd = `julia -t $t $(@__DIR__())/reco.jl`
  @info cmd
  run(cmd)
end

include("makeImages.jl")



