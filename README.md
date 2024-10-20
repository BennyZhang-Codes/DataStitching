# KomaHighOrder

This is a Julia toolbox for MR simulation that can incorporate dynamic field changes associated with the gradients throughout the sequence.
This is an extension of KomaMRI.jl (https://github.com/JuliaHealth/KomaMRI.jl), a Julia package for highly efficient MR simulations.
You may run the [demos](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo), to grab an idea of how this toolbox can be used to simulate MRI signals given a pulseq sequence and dynamic field changes measured using a field camera.

demos:

1. [SimRecon_SingleChannel](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/SimRecon_SingleChannel.jl): Simulation and reconstruction of a [fully-sampled single-shot spiral sequence (1 mm resolution)](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/Sequence/1mm_R1.seq) with [dynamic field changes](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/DynamicFields/1mm_R1.mat) and ΔB₀.
2. [Multi-echo_GRE](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/Muti-echo_GRE): estimating ΔB₀ map and coil-sensitivity map from the multi-echo GRE data in the [ISMRMRD](https://github.com/ismrmrd/ismrmrd) format. Additionally, the ME-GRE sequence can be modified within the [source code](https://github.com/BennyZhang-Codes/KomaHighOrder/blob/master/demo/Muti-echo_GRE/pulseq) ([Pulseq](https://github.com/pulseq/pulseq), MATLAB version).

If you use the toolbox, please consider citing the following abstract:

Zhang, Z., Auerbach, E. J., Bratch, A., Grant, A. N., Zhuo, Y., He, S., Chen, L., Ugurbil, K., Wu, X. "A stitching method for dynamic field monitoring using NMR probes", 2024 ISMRM, Singapore

### Copyright & License Notice

This software is copyrighted by the Regents of the University of Minnesota. It can be freely used for educational and research purposes by non-profit institutions and US government agencies only.
Other organizations are allowed to use this software only for evaluation purposes, and any further uses will require prior approval. The software may not be sold or redistributed without prior approval.
One may make copies of the software for their use provided that the copies, are not sold or distributed, are used under the same terms and conditions.
As unestablished research software, this code is provided on an "as is'' basis without warranty of any kind, either expressed or implied.
The downloading, or executing any part of this software constitutes an implicit agreement to these terms. These terms and conditions are subject to change at any time without prior notice.
