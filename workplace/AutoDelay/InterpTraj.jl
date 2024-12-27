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

plt_kspha_com(ksphaStitched, 
    InterpTrajTime(ksphaStitched, dt, dt*5, intermode=BSpline(Linear())), 
    dt)

kspha_akima  = InterpTrajTime(ksphaStitched, dt, dt*5, intermode=SteffenMonotonicInterpolation());
kspha_linear = InterpTrajTime(ksphaStitched, dt, dt*5, intermode=AkimaMonotonicInterpolation());
plt_kspha(kspha_linear - kspha_akima, dt)
    
kspha_akima  = InterpTrajTime(ksphaStitched, dt, dt*5, intermode=FritschButlandMonotonicInterpolation());
kspha_linear = InterpTrajTime(ksphaStitched, dt, dt*5, intermode=AkimaMonotonicInterpolation());
plt_kspha(kspha_linear - kspha_akima, dt)


"""
#=
    Interpolate using cubic Hermite splines. The breakpoints in arrays xbp and ybp are assumed to be sorted.
    Evaluate the function in all points of the array xeval.
    Methods:
        "Linear"                yuck
        "FiniteDifference"      classic cubic interpolation, no tension parameter
                                Finite difference can overshoot for non-monotonic data
        "Cardinal"              cubic cardinal splines, uses tension parameter which must be between [0,1]
                                cubin cardinal splines can overshoot for non-monotonic data
                                (increasing tension decreases overshoot)
        "Akima"                 monotonic - tangents are determined at each given point locally,
                                the curve obtained is close to a manually drawn curve, can overshoot for non-monotonic data
        "FritschCarlson"        monotonic - tangents are first initialized, then adjusted if they are not monotonic
                                can overshoot for non-monotonic data
        "FritschButland"        monotonic - faster algorithm (only requires one pass) but somewhat higher apparent "tension"
        "Steffen"               monotonic - also only one pass, results usually between FritschCarlson and FritschButland
    Sources:
        Akima (1970), "A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures", doi:10.1145/321607.321609
        Fritsch & Carlson (1980), "Monotone Piecewise Cubic Interpolation", doi:10.1137/0717021.
        Fritsch & Butland (1984), "A Method for Constructing Local Monotone Piecewise Cubic Interpolants", doi:10.1137/0905021.
        Steffen (1990), "A Simple Method for Monotonic Interpolation in One Dimension", http://adsabs.harvard.edu/abs/1990A%26A...239..443S

    Implementation based on http://bl.ocks.org/niclasmattsson/7bceb05fba6c71c78d507adae3d29417
=#

export
    LinearMonotonicInterpolation,
    FiniteDifferenceMonotonicInterpolation,
    CardinalMonotonicInterpolation,
    AkimaMonotonicInterpolation,
    FritschCarlsonMonotonicInterpolation,
    FritschButlandMonotonicInterpolation,
    SteffenMonotonicInterpolation
"""