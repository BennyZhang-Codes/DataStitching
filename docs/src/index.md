# HighOrderMRI.jl

This is a Julia toolbox for MR simulation and reconstruction that can incorporate dynamic field changes associated with the gradients throughout the sequence. This is an extension of [KomaMRI.jl](https://github.com/JuliaHealth/KomaMRI.jl) (a Julia package for highly efficient MR simulations) and [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl) (a Julia package for MRI reconstruction).

For MRI image reconstruction with field dynamics, we have built an extended signal encoding operator `HighOrderOp` to construct the signal equation. `HighOrderOp` inherits from `AbstractLinearOperator` in [LinearOperators.jl](https://github.com/JuliaSmoothOptimizers/LinearOperators.jl). Then, image reconstruction problem can be solved using algorithms from [RegularizedLeastSquares.jl](https://github.com/JuliaImageRecon/RegularizedLeastSquares.jl).

If you use the toolbox, please consider citing the following abstracts:

[1] Zhang, Z., Auerbach, E. J., Bratch, A., Grant, A. N., Zhuo, Y., He, S., Chen, L., Ugurbil, K., Wu, X. "A stitching method for dynamic field monitoring using NMR probes", 2024 ISMRM, Singapore

[2] Zhang, J., Zuo, Z., Xue, R., Zhuo, Y., Cushing, C., Bratch, A., Auerbach, E. J., Grant, A. N., Ugurbil, K., Wu, X., Zhang, Z. "A stitching method for dynamic field monitoring using NMR probes: validation in simulation and human experiments", 2025 ISMRM, Hawaii

### Features

* Support up to 2nd or 3rd order spherical harmonic terms.
* Support parrallel imaging and off-resonance correction with extended signal encoding operator `HighOrderOp`.
* Support the model-based synchronization delay estimation algorithm (Dubovan PI, Baron CA. 2023, [https://doi.org/10.1002/mrm.29460](https://doi.org/10.1002/mrm.29460)).
* GPU acceleration with `CUDA.jl` (only NVIDIA GPU has been tested). If the GPU memory is not enough, the calculation can be divided into blocks.

---

### Installation Guide

HighOrderMRI is compatible with Julia version 1.9.4. To get started with HighOrderMRI, users should first install Julia and consider using a code editor for a smoother coding experience.

To get the HighOrder package installed, execute the following Julia command:

```julia
import Pkg
Pkg.add(url="https://github.com/BennyZhang-Codes/HighOrderMRI.jl.git")
# or
Pkg.develop(url="https://github.com/BennyZhang-Codes/HighOrderMRI.jl.git")
```

---
