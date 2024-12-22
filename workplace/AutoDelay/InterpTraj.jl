using MAT
using Interpolations
using PyPlot


path = "$(@__DIR__)/workplace/AutoDelay/data"

grefile = "$(path)/syn_meas_MID00117_FID53005_pulseq_v0_gres6_1p0_standard.mat"
MRIfile = "$(path)/syn_meas_MID00115_FID53003_pulseq_v0_r4_1p0_standard.mat"
DFCfile = "$(path)/7T_1p0_200_r4.mat"


csm  = matread(grefile)["csm"];
mask = matread(grefile)["mask"];
b0   = matread(grefile)["b0"];

data       = matread(MRIfile)["data"];
matrixSize = matread(MRIfile)["matrixSize"];
FOV        = matread(MRIfile)["FOV"];

dt            = matread(DFCfile)["dt"];
ksphaStitched = matread(DFCfile)["ksphaStitched"];
ksphaStandard = matread(DFCfile)["ksphaStandard"];
delayStitched = matread(DFCfile)["delayStitched"];
delayStandard = matread(DFCfile)["delayStandard"];

bfieldStitched = matread(DFCfile)["bfieldStitched"];
bfieldStandard = matread(DFCfile)["bfieldStandard"];

plt_kspha(ksphaStitched, dt)
plt_bfield(bfieldStitched, dt)

plt_kspha_com(ksphaStitched, interp1_TrajTime(ksphaStitched, dt, dt*5), dt)



