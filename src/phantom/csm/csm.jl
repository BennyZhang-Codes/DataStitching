# This file is used to include all the necessary files for the Coil-Sensitivity Map (CSM) generating or loading.
include("Birdcage.jl")
include("Real_32cha.jl")

include("Fan.jl")
include("Rect.jl")
include("Ring.jl")
include("Gaussian.jl")
include("Gaussian_block.jl")

export csm_Birdcage
export csm_Real_32cha
export csm_Fan_binary
export csm_Rect_binary
export csm_Rect_gaussian
export csm_Ring_binary
export csm_Gaussian_grid
export csm_Gaussian_grid_block
export csm_Gaussian_grid_block_pha

include("load_csm.jl")
export load_csm