seq = read_seq("E:/skope/20240308/spiral_fov220_nx220_1.00.seq")
obj = brain_phantom2D(ss=10)
sys = Scanner()
sim_params = KomaMRICore.default_sim_params()
raw = HO_simulate(obj, seq, sys; sim_params=sim_params)


# Auxiliary function for reconstruction
function reconstruct_2d_image(raw::RawAcquisitionData)
    acqData = AcquisitionData(raw)
    acqData.traj[1].circular = false #Removing circular window
    C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
    acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
    Nx, Ny = raw.params["reconSize"][1:2]
    recParams = Dict{Symbol,Any}()
    recParams[:reconSize] = (Nx, Ny)
    recParams[:densityWeighting] = true
    rec = reconstruction(acqData, recParams)
    image3d  = reshape(rec.data, Nx, Ny, :)
    image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
    return image2d
end

# Perform reconstruction to get the image
image = reconstruct_2d_image(raw)
plot_image(image)