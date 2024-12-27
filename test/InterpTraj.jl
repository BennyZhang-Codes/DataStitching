using MAT
using Interpolations

nSample = 1000
nTerm   = 9
dt      = 1.    # dwell time for trajectory is 1

datatime = collect(0:2.5:nSample-1); # dwell time for data is 2.5
kspha  = rand(1000, 9);
kspha1 = InterpTrajTime(kspha, dt, dt*1);
kspha2 = InterpTrajTime(kspha, dt, 0., datatime);

plt_kspha(kspha, dt)
plt_kspha(kspha1, dt)
plt_kspha(kspha2, dt)
plt_kspha_com(kspha, kspha1, dt)


kspha_akima  = InterpTrajTime(kspha, dt, dt*1, intermode=SteffenMonotonicInterpolation());
kspha_linear = InterpTrajTime(kspha, dt, dt*1, intermode=AkimaMonotonicInterpolation());
plt_kspha(kspha_linear - kspha_akima, dt)
    
kspha_akima  = InterpTrajTime(kspha, dt, dt*1, intermode=FritschButlandMonotonicInterpolation());
kspha_linear = InterpTrajTime(kspha, dt, dt*1, intermode=AkimaMonotonicInterpolation());
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