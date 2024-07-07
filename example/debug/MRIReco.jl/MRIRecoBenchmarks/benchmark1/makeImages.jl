using Plots, Measures, HDF5, LinearAlgebra, RegularizedLeastSquares, DelimitedFiles, DataFrames
#pgfplotsx()


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
  header = ["Lib", "threads", "toeplitz", "oversampling", "r1", "r2", "r3", "r4"]
  data = readdlm(f_times, ',')
  df = DataFrame(data, vec(header))
  
  threads = unique(df[:,:threads])
  
  colors = [RGB(0.0,0.29,0.57), RGB(0.3,0.5,0.7), RGB(0.94,0.53,0.12), RGB(0.99,0.75,0.05)]
  
  Plots.scalefontsizes()
  Plots.scalefontsizes(0.99)
  
  ls = [:solid, :solid, :solid]
  shape = [:circle, :circle, :circle]
  lw = 3
  
  pl = Any[]
  
  for r in [1,4]
  
    p = plot(threads, df[df.Lib .== "BART", Symbol("r$r")],  #ylims=(0.0,1.2),
              lw=lw, xlabel = "# Threads", ylabel = "Time [s]", label="BART",
              legend = :topright,  
              title="R = $r", 
              shape=shape[1], ls=ls[1], 
              c=colors[1], msc=colors[1], mc=colors[1], ms=4, msw=2)
              
    plot!(p, threads, df[df.Lib .== "MRIReco" .&& df.toeplitz .== true, Symbol("r$r")], 
              label="MRIReco", lw=lw, shape=shape[2], ls=ls[2], 
              c=colors[2], msc=colors[2], mc=colors[2], ms=4, msw=2)
              
    plot!(p, threads, df[df.Lib .== "MRIReco" .&& df.toeplitz .== false, Symbol("r$r")], 
              label="MRIReco NonToeplitz", lw=lw, shape=shape[3], ls=ls[3], 
              c=colors[3], msc=colors[3], mc=colors[3], ms=4, msw=2)
              
    push!(pl, p)
  
  end
  
  p_ = plot(pl..., size=(900,400), layout=(1,2), bottom_margin=5mm, left_margin=5mm  )
             
  savefig(p_, @__DIR__() * "/reco/timings.svg")
 
  return p_
  
end


function makeImages()

  f_img  = @__DIR__() * "/reco/images.h5"
  imagesBart = reverse(abs.(h5read(f_img, "/recoBART")),dims=1)
  imagesMRIReco = reverse(abs.(h5read(f_img, "/recoMRIReco0")), dims=1)
  rf = size(imagesBart,3)
  
  
  for r = 1:rf
    imagesMRIReco[:,:,r] ./= maximum(imagesMRIReco[:,:,r])
    imagesBart[:,:,r] .= optimalScaling(imagesMRIReco[:,:,r], imagesBart[:,:,r])
  end
  m = 1
  
  plBART = Any[]
  for r = 1:rf
    push!(plBART, heatmap( imagesBart[:,:,r], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             title="R = $r", annotations = 
                (5,25, Plots.text((r==1) ? "BART" : "", :white, :left)) ) )
  end
  
  plBARTDiff = Any[]
  for r = 1:rf
    push!(plBARTDiff, heatmap( optimalScaling(imagesBart[:,:,1],imagesBart[:,:,r])  - imagesBart[:,:,1], 
             clim=(-0.1*m,0.1*m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             annotations = 
                (5,25, Plots.text((r==1) ? "Diff BART" : "", :white, :left)) ) )
  end
  
  plMRIReco = Any[]
  for r = 1:rf
    push!(plMRIReco, heatmap( imagesMRIReco[:,:,r], clim=(0,m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             annotations = 
                (5,25, Plots.text((r==1) ? "MRIReco.jl" : "", :white, :left)) ) )
  end
  
  plMRIRecoDiff = Any[]
  for r = 1:rf
    push!(plMRIRecoDiff, heatmap( optimalScaling(imagesMRIReco[:,:,1],imagesMRIReco[:,:,r]) 
                                  - imagesMRIReco[:,:,1], 
             clim=(-0.1*m,0.1*m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             annotations = 
                (5,25, Plots.text((r==1) ? "Diff MRIReco.jl" : "", :white, :left)) ) )
  end
  
  plDiff = Any[]
  for r = 1:rf
    push!(plDiff, heatmap( imagesMRIReco[:,:,r] - imagesBart[:,:,r], 
             clim=(-0.1*m,0.1*m), c=:viridis, 
             ticks=nothing, colorbar=nothing, 
             annotations = 
                (5,25, Plots.text((r==1) ? "MRIReco.jl-BART" : "", :white, :left)) ) )
  end
  
  p_ = plot(plBART..., plBARTDiff..., plMRIReco..., plMRIRecoDiff..., plDiff...,
             size=(900,900*5/4), layout=(5,4), left_margin = 0mm, right_margin=0mm )
 
  savefig(p_, @__DIR__() * "/reco/images.svg")
 
  return p_
end


makeImages()
makeTimings()





