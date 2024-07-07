# julia -t 4
ENV["OMP_NUM_THREADS"] = 4
using MRIReco, KomaHighOrder, MRICoilSensitivities

raw = RawAcquisitionData(ISMRMRDFile("E:/skope/KomaHighOrder/debug_recon/MRIReco.jl/xw_sp2d-1mm-r1_nominal.mrd"))

acqData = AcquisitionData(raw)
acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
Nx, Ny = raw.params["reconSize"][1:2]

encParams = MRIReco.getEncodingOperatorParams(;params...)
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, params)


T = Float32
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = T(1.e-2)
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
recParams[:correctionMap] =reshape(sensitivity, N, N, 1, Nc)
recParams[:encodingOps] = [fourierEncodingOp((N, N), tr[i], "explicit"; subsampleIdx=idx[i], correctionMap=reshape(sensitivity, N, N, 1, Nc),slice=1,encParams...) for i=1:numContr]

fourierEncodingOp((N, N), tr[i], "explicit"; subsampleIdx=idx[i], correctionMap=sensitivity[:,:])
rec = reconstruction(acqData, recParams)
image3d  = reshape(rec.data, Nx, Ny, :)
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
plot_image(image2d; title="iterative - CGNR", width=650, height=600)


import MRIReco.MRIOperators: ExplicitOp

shape = (Nx, Ny)
# sensitivity
acqDataCart = regrid(acqData, shape; cgnr_iter=3);
sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98);
# plot_image(abs.(rotr90(sensitivity[:,:,1,1])); title="sensitivity - espirit")


# encodingOps_simple
numContr = numContrasts(acqData)
tr = [trajectory(acqData,i) for i=1:numContr]
idx = acqData.subsampleIndices

F = fourierEncodingOp(shape, tr[1], "explicit", subsampleIdx=idx[1], correctionMap=sensitivity[:,:])

# ExplicitOp
nodes,times = MRIReco.kspaceNodes(tr[1]), MRIReco.readoutTimes(tr[1])
nrow = size(nodes,2)
ncol = prod(shape)
echoOffset = 0.0
correctionmap = sensitivity[:,:]
encOp = ExplicitOp{ComplexF64,Nothing,Function}(nrow, ncol, false, false
            , (res,x)->(res .= MRIReco.MRIOperators.produ(x, shape, Float64.(nodes), Float64.(times), echoOffset, ComplexF64.(correctionmap)))
            , nothing
            , (res,y)->(res .= MRIReco.MRIOperators.ctprodu(y, shape, Float64.(nodes), Float64.(times), echoOffset, ComplexF64.(correctionmap)))
            , 0,0,0, false, false, false, ComplexF64[], ComplexF64[])

T = Float32
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = T(1.e-2)
recParams[:iterations] = 30
recParams[:solver] = "cgnr"
recParams[:correctionMap] =sensitivity[:,:]
recParams[:encodingOps] = reshape([encOp], 1,1)

# fourierEncodingOp((N, N), tr[i], "explicit"; subsampleIdx=idx[i], correctionMap=sensitivity[:,:])
rec = reconstruction(acqData, recParams)
image3d  = reshape(rec.data, Nx, Ny, :);
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1];
plot_image(image2d; title="iterative - ExplicitOp_CGNR", width=650, height=600)






### NUFFT without density compensation
F = fourierEncodingOp(shape, tr[1], "fast", subsampleIdx=idx[1])
img = adjoint(F) *  kData(acqData,1,1,1,rep=1)
img = reshape(img, shape)
plot_image(abs.(img))