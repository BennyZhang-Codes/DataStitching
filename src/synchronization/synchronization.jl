include("InterpTraj.jl")
export InterpTrajTime

import DSP: conv               # using convolution to approximate derivative
include("FindDelay.jl")
include("FindDelay_v2.jl")
include("FindDelay_v2_multishot.jl")
export FindDelay, FindDelay_v2, FindDelay_v2_multishot


