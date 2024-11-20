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

csm_list = [:fan, :ring, :rect, :rect_gaussian, :birdcage, :real_32cha, :gaussian_grid, :gaussian_grid_block]
export csm_list

function load_csm(
    type::Symbol                    , 
    nX::Int64                       , 
    nY::Int64                       , 
    nCoil::Int64                    ;
    nRow                  = nothing ,
    nCol                  = nothing ,
    nBlock                = 3       ,
    overlap::Real         = 1       ,  # overlap between fan coils, for csm_Fan_binary
    relative_radius::Real = 1.5     ,  # relative radius of the coil, for csm_Birdcage
    verbose::Bool         = false
)
    @assert type in csm_list "type must be one of the following: $(csm_list)"
    if type == :fan
        csm = csm_Fan_binary(nX, nY, nCoil; overlap=overlap, verbose=verbose)
    elseif type == :ring
        csm = csm_Ring_binary(nX, nY, nCoil; overlap=overlap, verbose=verbose)
    elseif type == :rect
        csm = csm_Rect_binary(nX, nY, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif type == :rect_gaussian
        csm = csm_Rect_gaussian(nX, nY, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif type == :birdcage
        csm = csm_Birdcage(nX, nY, nCoil; relative_radius=Float64(relative_radius), verbose=verbose)
    elseif type == :real_32cha
        csm = csm_Real_32cha(nX, nY; verbose=verbose)
    elseif type == :gaussian_grid
        csm = csm_Gaussian_grid(nX, nY, nCoil; relative_radius=Float64(relative_radius), nRow=nRow, nCol=nCol, verbose=verbose)
    elseif type == :gaussian_grid_block
        csm = csm_Gaussian_grid_block(nX, nY, nCoil; relative_radius=Float64(relative_radius), nRow=nRow, nCol=nCol, nBlock=nBlock,verbose=verbose)
    end
    return csm
end