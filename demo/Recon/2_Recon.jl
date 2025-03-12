using MAT
using CUDA
CUDA.device!(0)

T             = Float64;
path          = joinpath(@__DIR__, "demo/Recon")
data_file     = "$(path)/data_1p0_200_r4.mat"

@info "data file: $(data_file)"
data          = matread(data_file);

csm           = data["gre_csm"];            # coil sensitivity map
b0            = data["gre_b0"];             # ΔB0 map
mask          = data["gre_mask"];           # mask

kdata         = data["kdata"];              # spiral k-space data
datatime      = data["datatime"];           # time stamps of k-space data
matrixSize    = data["matrixSize"];         # matrix size
FOV           = data["FOV"];                # field of view

k0_ecc        = data["k0_ecc"];             # b0 compensation of scanner from ECC model
dt            = data["dfc_dt"];             # dwell time of field dynamics
ksphaStitched = data["dfc_ksphaStitched"];  # coefficients of the field dynamics with data stitching method
ksphaStandard = data["dfc_ksphaStandard"];  # coefficients of the field dynamics with standard method
startStitched = data["dfc_startStitched"];  # start time of the field dynamics with data stitching method
startStandard = data["dfc_startStandard"];  # start time of the field dynamics with standard method

ksphaNominal  = data["ksphaNominal"];       # nominal kspace trajectory, kx, ky
startNominal  = data["startNominal"];       # start time of the nominal trajectory

tau_Nominal   = data["tau_Nominal"];        # synchronization delay between the nominal trajectory and the MRI data
tau_Stitched  = data["tau_Stitched"];       # synchronization delay between the stitched trajectory and the MRI data    
tau_Standard  = data["tau_Standard"];       # synchronization delay between the standard trajectory and the MRI data


# comparion of the field dynamics with or without synchronization
# plt_ksphas([InterpTrajTime(ksphaStitched, dt, startStitched, datatime), InterpTrajTime(ksphaStitched, dt, startStitched+tau_Stitched*dt, datatime)], dt)
# plt_ksphas([InterpTrajTime(ksphaStandard, dt, startStandard, datatime), InterpTrajTime(ksphaStandard, dt, startStandard+tau_Standard*dt, datatime)], dt)
# plt_ksphas([InterpTrajTime(ksphaNominal , dt, startNominal , datatime), InterpTrajTime(ksphaNominal , dt, startNominal +tau_Nominal *dt, datatime)], dt)

ksphaStitched = InterpTrajTime(ksphaStitched, dt, startStitched+tau_Stitched*dt, datatime);
ksphaStandard = InterpTrajTime(ksphaStandard, dt, startStandard+tau_Standard*dt, datatime);
ksphaNominal  = InterpTrajTime(ksphaNominal , dt, startNominal +tau_Nominal *dt, datatime);

# prepare some parameters for reconstruction
# 1. gridding 
Δx, Δy, Δz    = T.(FOV ./ matrixSize);
nX, nY, nZ    = matrixSize;
gridding      = Grid(nX, nY, nZ, Δx, Δy, Δz; exchange_xy=true, reverse_x=false, reverse_y=true)

# 2. sampling density 
weightStitched = SampleDensity(ksphaStitched'[2:3,:], (nX, nY));
weightStandard = SampleDensity(ksphaStandard'[2:3,:], (nX, nY));
weightNominal  = SampleDensity( ksphaNominal'[2:3,:], (nX, nY));

use_gpu = true;
verbose = false;
nBlock  = 20;

solver = "cgnr"; reg = "L2"; iter = 20; λ = 1e-9
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver


###########################################################################
# recon with field dynamics measured by our data stitching method
###########################################################################
kdata = data["kdata"];
kdata = kdata ./ exp.(1im.*k0_ecc)';

labelStitched = [    "Stitched",      "Stitched_wo_dB0"];
recons        = [         "111",                  "111"];
kdatas        = [         kdata,                  kdata];
weights       = [weightStitched,         weightStitched]; 
b0s           = [            b0,                  b0.*0];
imgStitched   = Array{Complex{T},3}(undef, nX, nY, length(labelStitched));
idxs           = collect(1:length(labelStitched));
for (idx, label, recon, kdata, weight, b0) in zip(idxs, labelStitched, recons, kdatas, weights, b0s)
    @info "[$(idx)] $(label)"
    HOOp = HighOrderOp(gridding, T.(ksphaStitched'), T.(datatime); recon_terms=recon, k_nominal=ksphaNominal[:, 2:4]', 
        nBlock=nBlock, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
    imgStitched[:, :, idx] = x;
    fig = plt_image(abs.(x); title=label, vmaxp=99.9)
    fig.savefig("$(path)/$(label).png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0.0)
end

###########################################################################
# recon with field dynamics measured by conventional single measurement
###########################################################################
kdata = data["kdata"];
kdata = kdata ./ exp.(1im.*k0_ecc)';

labelStandard = [    "Standard"];
recons        = [         "111"];
kdatas        = [         kdata];
weights       = [weightStandard]; 
b0s           = [            b0];
imgStandard   = Array{Complex{T},3}(undef, nX, nY, length(labelStandard));
idxs           = collect(1:length(labelStandard));
for (idx, label, recon, kdata, weight, b0) in zip(idxs, labelStandard, recons, kdatas, weights, b0s)
    @info "[$(idx)] $(label)"
    HOOp = HighOrderOp(gridding, T.(ksphaStandard'), T.(datatime); recon_terms=recon, k_nominal=ksphaNominal[:, 2:4]', 
        nBlock=nBlock, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
    imgStandard[:, :, idx] = x;
    fig = plt_image(abs.(x); title=label, vmaxp=99.9)
    fig.savefig("$(path)/$(label).png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0.0)
end


###########################################################################
# recon with nominal trajectory
###########################################################################
kdata = data["kdata"];

labelNominal  = [    "Nominal"];
recons        = [        "010"];
kdatas        = [        kdata];
weights       = [weightNominal]; 
b0s           = [           b0];
imgNominal    = Array{Complex{T},3}(undef, nX, nY, length(labelNominal));
idxs           = collect(1:length(labelNominal));
for (idx, label, recon, kdata, weight, b0) in zip(idxs, labelNominal, recons, kdatas, weights, b0s)
    @info "[$(idx)] $(label)"
    HOOp = HighOrderOp(gridding, T.(ksphaNominal'), T.(datatime); recon_terms=recon, k_nominal=ksphaNominal[:, 2:4]', 
        nBlock=nBlock, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    @time x = recon_HOOp(HOOp, Complex{T}.(kdata), Complex{T}.(weight), recParams);
    imgNominal[:, :, idx] = x;
    fig = plt_image(abs.(x); title=label, vmaxp=99.9)
    fig.savefig("$(path)/$(label).png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0.0)
end
