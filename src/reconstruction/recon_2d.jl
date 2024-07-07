function reconstruct_2d_image(raw::RawAcquisitionData; Nx=nothing, Ny=nothing)
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
    image3d  = reshape(rec.data, Nx, Ny, :)
    image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
    return image2d
end

function reconstruct_2d_image_multi_coil(raw::RawAcquisitionData; Nx=nothing, Ny=nothing)
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
    nX, nY, nZ, nCon, nCoil, nRep = size(rec)
    images = abs.(rec[:,:,1,1,:,1])
    return images
end

