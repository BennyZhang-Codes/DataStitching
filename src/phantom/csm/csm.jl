# This file is used to include all the necessary files for the Coil-Sensitivity Map (CSM) generating or loading.
include("Birdcage.jl")
include("Real_32cha.jl")

include("Fan.jl")
include("Rect.jl")
include("Ring.jl")

export csm_Birdcage
export csm_Real_32cha
export csm_Fan_binary
export csm_Rect_binary
export csm_Rect_gaussian
export csm_Ring_binary



function load_csm(
    type::Symbol                    , 
    nX::Int64                       , 
    nY::Int64                       , 
    nCoil::Int64                    ;
    nRow                  = nothing ,
    nCol                  = nothing ,
    overlap::Real         = 1       ,  # overlap between fan coils, for csm_Fan_binary
    relative_radius::Real = 1.5     ,  # relative radius of the coil, for csm_Birdcage
    verbose::Bool         = false
)
    @assert type in [:fan, :rect, :rect_gaussian, :birdcage, :real_32cha] "type must be one of the following: :fan, :rect, :rect_gaussian, :birdcage, :real_32cha"
    if type == :fan
        csm = csm_Fan_binary(nX, nY, nCoil; overlap=overlap, verbose=verbose)
    elseif type == :rect
        csm = csm_Rect_binary(nX, nY, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif type == :rect_gaussian
        csm = csm_Rect_gaussian(nX, nY, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif type == :birdcage
        csm = csm_Birdcage(nX, nY, nCoil; relative_radius=Float64(relative_radius), verbose=verbose)
    elseif type == :real_32cha
        csm = csm_Real_32cha(nX, nY; verbose=verbose)
    end
    return csm
end