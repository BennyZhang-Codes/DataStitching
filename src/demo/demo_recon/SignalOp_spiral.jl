# julia -t 4
device!(1)
name.(devices())

ENV["OMP_NUM_THREADS"] = 4
using MRIReco, KomaHighOrder, MRICoilSensitivities, ProgressMeter

raw = demo_raw("000")

Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false
shape = (Nx, Ny)

HOOp = SignalOp(shape, acqData.traj[1]; Nblocks=3)

T = Float32
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = T(1.e-2)
recParams[:iterations] = 1
recParams[:solver] = "cgnr"
recParams[:encodingOps] = reshape([HOOp], 1,1)

@time rec = reconstruction(acqData, recParams);
image3d  = reshape(rec.data, Nx, Ny, :);
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1];
plot_image(image2d; title="iterative_CG-SignalOp", width=650, height=600)



