# julia -t 4
# using CUDA
# device!(1) 

using MRIReco, KomaHighOrder, MRICoilSensitivities, PlotlyJS, MAT
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral"
BHO_simu = "000"
raw = demo_raw(BHO_simu)
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq()
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1)

tr_skope = Trajectory(K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);

#######################################################################################
# HighOrderOp 
# [1] Simu: 111, Reco: 000
# [2] Simu: 111, Reco: 100
# [3] Simu: 111, Reco: 010
# [4] Simu: 111, Reco: 001
# [5] Simu: 111, Reco: 011
# [6] Simu: 111, Reco: 101
# [7] Simu: 111, Reco: 110
# [8] Simu: 111, Reco: 111
#######################################################################################
BHO_simu = "111";raw = demo_raw(BHO_simu);
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);acqData.traj[1].circular = false;shape = (Nx, Ny);

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
BHO_recos = ["000", "100", "010", "001", "011", "101", "110", "111"]
imgs = Array{ComplexF32,3}(undef, Nx, Ny, length(BHO_recos))
for idx in eachindex(BHO_recos)
    @info "Simu: $(BHO_simu), Reco: $(BHO_recos[idx])"
    BHO_reco = BHO_recos[idx]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    imgs[:,:, idx] = rec.data[:,:]
end

imgs_111 = abs.(imgs)
imgs_111_error = Array{Float32,3}(undef, size(imgs_111))
for idx in eachindex(BHO_recos)
    imgs_111_error[:,:, idx] = imgs_111[:,:, idx] - imgs_111[:,:, end]
end

subplot_titles = ["Reco: $t" for t in BHO_recos]
title="HighOrderOp, Simu: $(BHO_simu)"

p_111       = plot_imgs(imgs_111, subplot_titles; title=title, width=1300, height=200)
p_111_error = plot_imgs(imgs_111_error, subplot_titles; title=title*", error map", width=1300, height=200)

savefig(p_111,       dir*"/HighOrderOp_Simu_111.svg", width=1300, height=200,format="svg")
savefig(p_111_error, dir*"/HighOrderOp_Simu_111_errormap.svg", width=1300, height=200,format="svg")
MAT.matwrite(dir*"/HighOrderOp_Simu_111.mat", Dict("imgs"=>imgs_111, "imgs_error"=>imgs_111_error, "BHO"=>BHO_recos))
#######################################################################################
# HighOrderOp 
# [1] Simu: 000, Reco: 000
# [2] Simu: 000, Reco: 100
# [3] Simu: 000, Reco: 010
# [4] Simu: 000, Reco: 001
# [5] Simu: 000, Reco: 011
# [6] Simu: 000, Reco: 101
# [7] Simu: 000, Reco: 110
# [8] Simu: 000, Reco: 111
#######################################################################################
BHO_simu = "000";raw = demo_raw(BHO_simu);
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);acqData.traj[1].circular = false;shape = (Nx, Ny);

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
BHO_recos = ["000", "100", "010", "001", "011", "101", "110", "111"]
imgs = Array{ComplexF32,3}(undef, Nx, Ny, length(BHO_recos))
for idx in eachindex(BHO_recos)
    @info "Simu: $(BHO_simu), Reco: $(BHO_recos[idx])"
    BHO_reco = BHO_recos[idx]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    imgs[:,:, idx] = rec.data[:,:]
end

imgs_000 = abs.(imgs)
imgs_000_error = Array{Float32,3}(undef, size(imgs_000))
for idx in eachindex(BHO_recos)
    imgs_000_error[:,:, idx] = imgs_000[:,:, idx] - imgs_000[:,:, 1]
end

subplot_titles = ["Reco: $t" for t in BHO_recos]
title="HighOrderOp, Simu: $(BHO_simu)"

p_000       = plot_imgs(imgs_000, subplot_titles; title=title, width=1300, height=200)
p_000_error = plot_imgs(imgs_000_error, subplot_titles; title=title*", error map", width=1300, height=200)

savefig(p_000,       dir*"/HighOrderOp_Simu_000.svg", width=1300, height=200,format="svg")
savefig(p_000_error, dir*"/HighOrderOp_Simu_000_errormap.svg", width=1300, height=200,format="svg")
MAT.matwrite(dir*"/HighOrderOp_Simu_000.mat", Dict("imgs"=>imgs_000, "imgs_error"=>imgs_000_error, "BHO"=>BHO_recos))







