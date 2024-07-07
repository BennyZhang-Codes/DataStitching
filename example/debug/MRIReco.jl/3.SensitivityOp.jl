using HDF5
using MRIReco
using KomaHighOrder, PlotlyJS

path = @__DIR__
path = path*"/src/demo/demo_recon/radial-ISMRM_reproducibility_challenge_1"
f_sensitivity  = path*"/data/sensitivitiesMRIReco.h5"
filename = path*"/data/rawdata_brain_radial_96proj_12ch.h5";

data = permutedims(h5read(filename, "rawdata"),[3,2,1,4]);
traj = permutedims(h5read(filename, "trajectory"),[3,2,1]);
N = 300;
Nc = 12;
T = Float32;
toeplitz = parse(Int,get(ENV,"TOEPLITZ","0"));
oversamplingFactor = parse(Float64,get(ENV,"OVERSAMPLING","1.25"));

tr = Trajectory(T.(reshape(traj[1:2,:,:],2,:) ./ N), 96, 512, circular=false);
dat = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1] = 1.e8.*reshape(data,:,12);
acqData = AcquisitionData(tr, dat, encodingSize=(N,N));

sensitivity = h5read(f_sensitivity, "/sensitivity");

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 20
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = toeplitz == 1
params[:oversamplingFactor] = oversamplingFactor
params[:senseMaps] = Complex{T}.(reshape(sensitivity, N, N, 1, Nc))

@info "undersampled reco"
rf = [1,2,3,4]
img_cg = Array{ComplexF32,3}(undef,N,N,4)
for i = eachindex(rf)
  @info "r=$(rf[i])"
  global acqDataSub = convertUndersampledData(sample_kspace(acqData, T(rf[i]), "regular"))
  img_cg[:,:,i] = reconstruction(acqDataSub, params).data
end

img = ComplexF64.(img_cg[:,:,1])
plot_image(abs.(img))

params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = encodingSize(acqData)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 20
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = toeplitz == 1
params[:oversamplingFactor] = oversamplingFactor
params[:senseMaps] = Complex{T}.(reshape(sensitivity, N, N, 1, Nc));

numContr, numChan = numContrasts(acqData), numChannels(acqData)
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, params)



smaps = senseMaps[:,:,1,:]
smaps = reshape(ComplexF64.(smaps),:,numChan)
S = SensitivityOp(smaps,1)

numVox, numChan = size(smaps)

x = vec(img)

numContr = 1
x_ = reshape(x,numVox,numContr)
y_ = Array{ComplexF64,3}(undef,numVox,numChan,numContr)

@assert size(smaps) == (size(y_,1), size(y_,2))

@time @inbounds for i ∈ CartesianIndices(y_)
    y_[i] = x_[i[1],i[3]] * smaps[i[1],i[2]]
  end


plot_image(abs.(reshape(y_, N, N, numChan)[:,:,1]))

