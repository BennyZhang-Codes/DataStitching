# julia -t 4
# using CUDA
# device!(1) 

using KomaHighOrder
using MRIReco, MRICoilSensitivities, PlotlyJS, MAT
simtype = SimType(B0=false, T2=true, ss=5)
BHO = BlochHighOrder("111")
skope_method = "Stitched"   # :Stitched or :Standard
dir = "$(@__DIR__)/src/demo/demo_recon/HighOrderOp_spiral/results_$skope_method/$(simtype.name)"; if ispath(dir) == false mkdir(dir) end


raw = demo_raw(BHO; simtype=simtype)
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;
shape = (Nx, Ny);

hoseq = demo_hoseq(skope_method=Symbol(skope_method))
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1)

tr_skope = Trajectory(K_skope_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);
tr_nominal = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);
# plot_traj2d(Trajectory(K_skope_adc'[2:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false))
# plot_traj2d(Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false))
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
BHO = BlochHighOrder("111"); raw = demo_raw(BHO; simtype=simtype)
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
imgs = Array{ComplexF32,3}(undef, Nx, Ny, length(BHO_recos));
for idx in eachindex(BHO_recos)
    @info "Simu: $(BHO.name), Reco: $(BHO_recos[idx])"
    BHO_reco = BHO_recos[idx]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    imgs[:,:, idx] = rec.data[:,:]
end

imgs_111 = abs.(imgs);
imgs_111_error = Array{Float32,3}(undef, size(imgs_111));
for idx in eachindex(BHO_recos)
    imgs_111_error[:,:, idx] = imgs_111[:,:, idx] - imgs_111[:,:, end]
end

subplot_titles = ["Reco: $t" for t in BHO_recos]
title="HighOrderOp, Simu: $(BHO.name)"

width = 1200 
height = 160
p_111       = plot_imgs(imgs_111, subplot_titles; title=title, width=width+100, height=height+40)
p_111_error = plot_imgs(imgs_111_error, subplot_titles; title=title*", error map", width=width+100, height=height+40)

savefig(p_111,       dir*"/HighOrderOp_Simu_$(BHO.name).svg", width=width+100, height=height+40,format="svg")
savefig(p_111_error, dir*"/HighOrderOp_Simu_$(BHO.name)_errormap.svg", width=width+100, height=height+40,format="svg")
MAT.matwrite(dir*"/HighOrderOp_Simu_$(BHO.name).mat", Dict("imgs"=>imgs_111, "imgs_error"=>imgs_111_error, "BHO"=>BHO_recos))
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
BHO = BlochHighOrder("000"); raw = demo_raw(BHO; simtype=simtype)
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
imgs = Array{ComplexF32,3}(undef, Nx, Ny, length(BHO_recos));
for idx in eachindex(BHO_recos)
    @info "Simu: $(BHO_simu.name), Reco: $(BHO_recos[idx])"
    BHO_reco = BHO_recos[idx]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    imgs[:,:, idx] = rec.data[:,:]
end

imgs_000 = abs.(imgs);
imgs_000_error = Array{Float32,3}(undef, size(imgs_000));
for idx in eachindex(BHO_recos)
    imgs_000_error[:,:, idx] = imgs_000[:,:, idx] - imgs_000[:,:, 1]
end

subplot_titles = ["Reco: $t" for t in BHO_recos]
title="HighOrderOp, Simu: $(BHO_simu.name)"

p_000       = plot_imgs(imgs_000, subplot_titles; title=title, width=1300, height=200)
p_000_error = plot_imgs(imgs_000_error, subplot_titles; title=title*", error map", width=1300, height=200)

savefig(p_000,       dir*"/HighOrderOp_Simu_$(BHO_simu.name).svg", width=1300, height=200,format="svg")
savefig(p_000_error, dir*"/HighOrderOp_Simu_$(BHO_simu.name)_errormap.svg", width=1300, height=200,format="svg")
MAT.matwrite(dir*"/HighOrderOp_Simu_$(BHO_simu.name).mat", Dict("imgs"=>imgs_000, "imgs_error"=>imgs_000_error, "BHO"=>BHO_recos))







