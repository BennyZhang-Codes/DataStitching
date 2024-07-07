function recon_2d(raw::RawAcquisitionData; Nx=nothing, Ny=nothing, outtype::Symbol=:mag, verbose=false)
    @assert outtype in [:mag, :pha, :raw] "outtype must be :mag, :pha, :raw"
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
    if isnothing(Nx) || isnothing(Ny)
        Nx, Ny = raw.params["reconSize"][1:2]
    end
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    # nX, nY, nZ, nCon, nCoil, nRep = size(rec)
    # size_info = Dict("nX"=>nX, "nY"=>nY, "nZ"=>nZ, "nCon"=>nCon, "nCoil"=>nCoil, "nRep"=>nRep) 
    dims = ["nX", "nY", "nZ", "nCon", "nCoil", "nRep"]

    dims_squeezed = dims[findall(size(rec).>1)]
    size_squeezed = size(rec)[findall(size(rec).>1)]

    info_string = "\n  size = [ "
    for idx in eachindex(dims_squeezed)
        if idx != 1
            info_string *= " | "
        end
        info_string *= "$(dims_squeezed[idx])=$(size_squeezed[idx])"
    end
    info_string *= " ]"
    if verbose
        @info "recon_2d: $(info_string)\n  type = $(String(outtype))"
    end

    if outtype == :mag
        images = abs.(rec)
    elseif outtype == :pha
        images = angle.(rec)
    else
        images = rec
    end

    for dim in reverse(sort(findall(size(rec).==1)))
        images = dropdims(images, dims=dim)
    end
    return images
end