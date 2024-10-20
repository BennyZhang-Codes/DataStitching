# MRIRecoBenchmarks
Benchmarks for MRI Reconstruction Frameworks

## Installation 

In order to run the benchmark suite you need to install Julia (version 1.8 or higher) and BART (version 0.7 or higher). The you need to udjust the BART path in the file `configuration.jl`, e.g. we set it to
```julia
ENV["TOOLBOX_PATH"] = "/opt/software/bart-0.7.00"
```
Next you need to install various Julia packages. To do so, go to the root folder off this repository and call Julia with
```
julia --project=. 
```
Then run the command 
```
julia> using Pkg; Pkg.instantiate()
```
which will install the required packages. Finally you need to install an unregistered Julia package manually by running
```
julia> Pkg.add(url="https://github.com/JakobAsslaender/BartIO.jl")
```
Now, everything is in place to run the benchmark suite. Next time you startup Julia, you just need `julia --project=.` since all packages are already installed.

## Run

To run the bachmark suite call
```
julia> include("runBenchmarks.jl")
```
Depending on your hardware, this can last some time (more than 10 minutes). The benchmarks will generate various images in the subfolders (e.g. benchmark1/reco/timings.svg). These image are included in the `Readme.md` files that you can watch by clicking on one subfolder on the GitHub web page.

## Hardware

The benchmark results shown on the webpage were performed on a computer with 1024 GB of main memory and an AMD EPYC 7702 CPU.
