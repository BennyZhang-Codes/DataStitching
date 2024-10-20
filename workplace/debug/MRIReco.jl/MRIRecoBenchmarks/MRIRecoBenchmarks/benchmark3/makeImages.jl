using Plots, Measures, HDF5, LinearAlgebra, RegularizedLeastSquares, DelimitedFiles, DataFrames


function optimalScaling(I,Ireco)
  N = length(I)

  # This is a little trick. We usually are not interested in simple scalings
  # and therefore "calibrate" them away
  alpha = norm(Ireco)>0 ? (dot(vec(I),vec(Ireco))+dot(vec(Ireco),vec(I))) /
          (2*dot(vec(Ireco),vec(Ireco))) : 1.0
  I2 = Ireco.*alpha

  return I2
end



function makeTimings()
  f_times = @__DIR__() * "/reco/recoTimes.csv"
  header = ["Lib", "threads", "time"]
  data = readdlm(f_times, ',')
  df = DataFrame(data, vec(header))
  
  threads = unique(df[:,:threads])
  
  colors = [RGB(0.0,0.29,0.57), RGB(0.3,0.5,0.7), RGB(0.94,0.53,0.12), RGB(0.99,0.75,0.05)]
  
  Plots.scalefontsizes()
  Plots.scalefontsizes(0.99)
  
  ls = [:solid, :solid, :solid]
  shape = [:circle, :circle, :circle]
  lw = 3
  

  p = plot(threads, df[df.Lib .== "BART", :time],  #ylims=(0.0,1.2),
              lw=lw, xlabel = "# Threads", ylabel = "Time [s]", label="BART",
              legend = :topright,  #yaxis = :log,
              shape=shape[1], ls=ls[1], 
              c=colors[1], msc=colors[1], mc=colors[1], ms=4, msw=2,
              size=(900,400), bottom_margin=5mm, left_margin=5mm )
              
  plot!(p, threads, df[df.Lib .== "MRIReco", :time], 
              label="MRIReco", lw=lw, shape=shape[2], ls=ls[2], 
              c=colors[2], msc=colors[2], mc=colors[2], ms=4, msw=2)
              
             
  savefig(p, @__DIR__() * "/reco/timings.svg")
 
  return p
  
end


function makeImages()

  f_img  = @__DIR__() * "/reco/images.h5"
  
  sensitivity = reverse(abs.(h5read(f_img, "/sensitivity")),dims=1)
  sensitivityBART = reverse(abs.(h5read(f_img, "/sensitivityBART")),dims=1)
  sensitivityMRIReco = reverse(abs.(h5read(f_img, "/sensitivityMRIReco")), dims=1)

  sensitivity./= maximum(sensitivity)
  sensitivityBART .= optimalScaling(sensitivity, sensitivityBART)
  sensitivityMRIReco .= optimalScaling(sensitivity, sensitivityMRIReco)

  m = 1
  slice = 64
  R = size(sensitivity,4)

  plTruth = Any[]
  for r = 1:R
    push!(plTruth, heatmap( sensitivity[:,:,slice,r], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             title="coil $r", annotations = 
                (5,25, Plots.text((r==1) ? "Truth" : "", :white, :left)) ) )
  end

  plBART = Any[]
  for r = 1:R
    push!(plBART, heatmap( sensitivityBART[:,:,slice,r], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             annotations = 
                (5,25, Plots.text((r==1) ? "BART" : "", :white, :left)) ) )
  end

  plMRIReco = Any[]
  for r = 1:R
    push!(plMRIReco, heatmap( sensitivityMRIReco[:,:,slice,r,1], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             annotations = 
                (5,25, Plots.text((r==1) ? "MRIReco" : "", :white, :left)) ) )
  end
  
                
  p_ = plot(plTruth..., plBART..., plMRIReco...,
             size=(900,400), layout=(3,R), left_margin = 0mm, right_margin=0mm )
 
  savefig(p_, @__DIR__() * "/reco/images.svg")
 
  return p_
end


makeImages()
makeTimings()





