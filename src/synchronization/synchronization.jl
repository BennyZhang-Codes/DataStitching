include("InterpTraj.jl")
export InterpTrajTime

import DSP: conv               # using convolution to approximate derivative
include("FindDelay.jl")
include("FindDelay_v2.jl")
export FindDelay, FindDelay_v2


