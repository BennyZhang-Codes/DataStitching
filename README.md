# KomaHighOrder.jl

This is a Julia toolbox for MR simulation that can incorporate dynamic field changes associated with the gradients throughout the sequence.
This is an extension of KomaMRI.jl (https://github.com/JuliaHealth/KomaMRI.jl), a Julia package for highly efficient MR simulations.
This repository provides a set of demo scripts that allow users to perform simulations and reconstructions using both our data stitching method and the standard approach.
You may run the [demos](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo), to grab an idea of how this toolbox can be used to simulate MRI signals given a pulseq sequence and dynamic field changes measured using a field camera.

If you use the toolbox, please consider citing the following abstract:

Zhang, Z., Auerbach, E. J., Bratch, A., Grant, A. N., Zhuo, Y., He, S., Chen, L., Ugurbil, K., Wu, X. "A stitching method for dynamic field monitoring using NMR probes", 2024 ISMRM, Singapore

## Methods

To enable the sequence simulation with higher spatial order terms of time-resolved field dynamics, we extended **KomaMRI.jl**  by including the fields up to second order in the calculation of $B_{z}$:

$$
B_{z}(t)=\sum_{l}b_l(t)h_l(\boldsymbol{r}) + \frac{\Delta\omega(\boldsymbol{r})}{\gamma},
$$

where $b_l(t)$ are the fields derived from the field dynamics for spherical harmonic basis functions $h_l(\boldsymbol{r})$. $l$ is the index denotes different spherical harmonic terms, $\Delta\omega(\boldsymbol{r})$ is the off-resonance field, $\gamma$ is the gyromagnetic ratio, and $\boldsymbol{r}$ is the position vector.
For image reconstruction, we implemented an extended signal encoding model based on **MRIReco.jl**. The measured or simulated MRI signal $s$ received in the coil $p$ at time $t$ can be describe by the model as follows:

$$
{s}_p(t)={\sum_n}{c}_p\left({\boldsymbol r}_n\right)m\left({\boldsymbol r}_n\right) {e}^{i \phi}{e}^{i\Delta \omega \left({\boldsymbol r}_n\right){t}}
$$

$$
\phi = {\sum_{l}}{k}_l\left({t} \right){h}_l\left({\boldsymbol r}_n\right)
$$

where $k_l(t)$ are the coefficients measured by NMR probes and $m$ represents magnetization. Using this model, images can be reconstructed with measured field dynamics and static off-resonance by iterative SENSE algorithm.

## Demo

1. [Sim and Recon for single-channel](https://github.com/BennyZhang-Codes/KomaHighOrder/tree/master/demo/SingleChannel): Simulation and reconstruction of a fully-sampled single-shot spiral sequence (1 mm resolution) [[.seq file]](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/SingleChannel/1mm_R1.seq) with field dynamics and ΔB₀. Reconstruction is based on a extended signal encoding model, which includes the field dynamics (up to second-order) and off-resonance.
2. [Sim and Recon for multi-channel](https://github.com/BennyZhang-Codes/KomaHighOrder/tree/master/demo/MultiChannel): Simulation and reconstruction of a under-sampled single-shot spiral sequence (1 mm resolution, R=4) [[.seq file]](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/MultiChannel/xw_sp2d_7T-1mm-200-r4-noSync-fa90.seq) with field dynamics and ΔB₀.
3. [Multi-echo Gradient Echo (ME-GRE)](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/Muti-echo_GRE): estimating ΔB₀ map (MRIFieldmaps.jl) and coil-sensitivity map (ESPIRiT) from the ME-GRE data in the [ISMRMRD](https://github.com/ismrmrd/ismrmrd) format. Additionally, the ME-GRE sequence can be modified within the [source code](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/Muti-echo_GRE/pulseq) ([Pulseq](https://github.com/pulseq/pulseq), MATLAB version).

### Required dependencies:

The current version is mainly based on two other packages: KomaMRI.jl (version: 0.8.0) and MRIReco.jl (version: 0.8.0).

- `Julia` 1.9.4
- `KomaMRI` 0.8.0
- `MRIReco` 0.8.0
- `RegularizedLeastSquares` 0.10.0
- `MRIFieldmaps` 0.0.3
- `PyPlot` 2.11.5

### Copyright & License Notice

This software is copyrighted by the Regents of the University of Minnesota. It can be freely used for educational and research purposes by non-profit institutions and US government agencies only.
Other organizations are allowed to use this software only for evaluation purposes, and any further uses will require prior approval. The software may not be sold or redistributed without prior approval.
One may make copies of the software for their use provided that the copies, are not sold or distributed, are used under the same terms and conditions.
As unestablished research software, this code is provided on an "as is'' basis without warranty of any kind, either expressed or implied.
The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.
