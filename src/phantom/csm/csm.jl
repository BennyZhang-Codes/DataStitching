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
    csmtype::Symbol, 
    Nx::Int64, 
    Ny::Int64, 
    nCoil::Int64;
    nRow            = nothing,
    nCol            = nothing,
    overlap::Real        =1,       # overlap between fan coils, for csm_Fan_binary
    relative_radius::Real=1.5,     # relative radius of the coil, for csm_Birdcage
    verbose::Bool=false
)
    @assert csmtype in [:fan, :rect, :birdcage, :real_32cha] "csmtype must be one of the following: :fan, :rect, :birdcage, :real_32cha"
    if csmtype == :fan
        csm = csm_Fan_binary(Nx, Ny, nCoil; overlap=overlap, verbose=verbose)
    elseif csmtype == :rect
        csm = csm_Rect_binary(Nx, Ny, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif csmtype == :rect_gaussian
        csm = csm_Rect_gaussian(Nx, Ny, nCoil; verbose=verbose, nRow=nRow, nCol=nCol)
    elseif csmtype == :birdcage
        csm = csm_Birdcage(Nx, Ny, nCoil; relative_radius=Float64(relative_radius), verbose=verbose)
    elseif csmtype == :real_32cha
        csm = csm_Real_32cha(Nx, Ny; verbose=verbose)
    end
    return csm
end