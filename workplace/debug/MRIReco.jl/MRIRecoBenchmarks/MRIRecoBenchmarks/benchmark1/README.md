## Benchmark 1 -  SENSE Reconstruction

The first benchmark runs the [ISMRM Reproducibility Challenges](https://ismrm.github.io/rrsg/) considering the brain data. A complete MRIReco.jl based implementation can be found [here](https://github.com/MagneticResonanceImaging/ISMRM_RRSG).

The reconstruction algorithm being run is the iterative SENSE reconstruction based on [this](https://doi.org/10.1002/mrm.1241) paper. It uses a conjugative gradient algorithm and appropriate gridding operators. The results are run for reduction factors R=1-4.

The MRIReco.jl implementation runs the benchmark twice. Once with the Toeplitz optimization for the normal matrix. A second time without that optimization but with ordinary NFFT based gridding operator using accuracy parameters tuned for maximum performance.

## Reconstruction Results

An accuracy comparison is provided in the following figure.

![Reconstruction Results](reco/images.svg?raw=true "Reconstruction Results")

## Benchmark Results

Benchmark results are given in the following figure.

![Benchmark Results](reco/timings.svg?raw=true "Benchmark Results")
