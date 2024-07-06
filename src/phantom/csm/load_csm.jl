


function load_csm(csmtype::Symbol, Nx::Int64, Ny::Int64, nCoils::Int64; verbose::Bool=false)
    @assert csmtype in [:fan, :rect, :birdcage, :real] "csmtype must be one of the following: :fan, :rect, :birdcage, :real_32cha"
    if csmtype == :fan
        csm = csm_Fan_binary(Nx, Ny, nCoils; overlap=overlap, verbose=verbose)
    elseif csmtype == :rect
        csm = csm_Rect_binary(Nx, Ny, nCoils; verbose=verbose)
    elseif csmtype == :birdcage
        csm = csm_Birdcage_binary(Nx, Ny, nCoils; relative_radius=Float64(relative_radius), verbose=verbose)
    elseif csmtype == :real_32cha
        csm = csm_Real_32cha(Nx, Ny; verbose=verbose)
    end
    return csm
end