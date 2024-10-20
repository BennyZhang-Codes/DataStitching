using MRIReco, KomaHighOrder, MRICoilSensitivities

# read raw data(*.mrd) and recon
raw = RawAcquisitionData(ISMRMRDFile("E:/skope/KomaHighOrder/debug_recon/MRIReco.jl/xw_sp2d-1mm-r1_nominal.mrd"))
# raw = RawAcquisitionData(ISMRMRDFile("E:/skope/KomaHighOrder/debug_recon/MRIReco.jl/meas_MID00569_FID26776_PM_0_5iso_tb8_p2_d0.mrd"))

acqData = AcquisitionData(raw)
acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
Nx, Ny = raw.params["reconSize"][1:2]

##### NUFFT
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
rec = reconstruction(acqData, recParams)
image3d  = reshape(rec.data, Nx, Ny, :)
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
plot_image(image2d; title="NUFFT", width=650, height=600)

##### iterative - CGNR
T = Float32
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = T(1.e-2)
recParams[:iterations] = 100
recParams[:solver] = "cgnr"

rec = reconstruction(acqData, recParams)
image3d  = reshape(rec.data, Nx, Ny, :)
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
plot_image(image2d; title="iterative - CGNR", width=650, height=600)


##### SENSE (1 coil)
raw = RawAcquisitionData(ISMRMRDFile("E:/skope/KomaHighOrder/debug_recon/MRIReco.jl/xw_sp2d-1mm-r1_nominal.mrd"))
acqData = AcquisitionData(raw)
acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C

# s = acqData.kdata[1,1,1]
# acqData.kdata[1,1,1] = [s s]
N = 149; Nc = 1
acqDataCart = regrid(acqData, (N,N); cgnr_iter=3)
sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98)
# s = [sensitivity[:,:,1,i] for i=1:2]
# s = [s[1] s[2]]
# plot_image(abs.(rotr90(s)); title="sensitivity - espirit")
plot_image(abs.(rotr90(sensitivity[:,:,1,1])); title="sensitivity - espirit")
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 100
params[:solver] = "cgnr"
params[:toeplitz] = false
params[:oversamplingFactor] = 1.25
params[:senseMaps] = Complex{T}.(reshape(sensitivity, N, N, 1, Nc))

# params[:iterations] = 20
params[:relTol] = 0.0
img_cg = reconstruction(acqData, params).data
plot_image(abs.(img_cg[:,:]); title="R = 1")


