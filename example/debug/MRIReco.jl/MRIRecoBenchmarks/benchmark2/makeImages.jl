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
              legend = :topright,  
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
  
  imFully = reverse(abs.(h5read(f_img, "/imFully")),dims=1)
  imagesBart = reverse(abs.(h5read(f_img, "/recoBART")),dims=1)
  imagesMRIReco = reverse(abs.(h5read(f_img, "/recoMRIReco")), dims=1)

  imFully[:,:,:] ./= maximum(imFully[:,:,:])
  imagesBart[:,:,:] .= optimalScaling(imFully[:,:,:], imagesBart[:,:,:])
  imagesMRIReco[:,:,:] .= optimalScaling(imFully[:,:,:], imagesMRIReco[:,:,:])

  m = 1
  slice = 80
  
  plFully = heatmap( imFully[:,:,slice], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             title="FullySampled" )   
  
  plBART = heatmap( imagesBart[:,:,slice], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             title="BART" ) 
             
  plMRIReco = heatmap( imagesMRIReco[:,:,slice], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             title="MRIReco.jl" ) 
                
                
  p_ = plot(plFully, plBART, plMRIReco,
             size=(900,300), layout=(1,3), left_margin = 0mm, right_margin=0mm )
 
  savefig(p_, @__DIR__() * "/reco/images.svg")
 
  return p_
end


makeImages()
makeTimings()





