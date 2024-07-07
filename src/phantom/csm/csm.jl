# This file is used to include all the necessary files for the Coil-Sensitivity Map (CSM) generating or loading.
include("Birdcage.jl")
include("binary.jl")
include("Real_32cha.jl")

export csm_Birdcage
export csm_Real_32cha
export csm_Fan_binary
export csm_Rect_binary

include("load_csm.jl")
export load_csm



