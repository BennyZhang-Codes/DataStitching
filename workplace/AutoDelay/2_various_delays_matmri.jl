using KomaHighOrder
using MAT
using PyPlot
using MRIReco
using CUDA

CUDA.device!(2)

T=Float64;
path = "$(@__DIR__)/workplace/AutoDelay/comparison_matmri"

data_file = "$(path)/data_images.mat"
traj_file = "$(path)/data_traj_spiral_R4.mat"
data_mat = matread(data_file)
traj_mat = matread(traj_file)

nSli = 1;
x    = data_mat["X"][:,:,nSli];
y    = data_mat["Y"][:,:,nSli];
z    = data_mat["Z"][:,:,nSli];
csm  = data_mat["Crcvr"][:,:,nSli,:];
b0   = data_mat["b0map"][:,:,nSli] ./ (2π);
im0  = data_mat["im0"][:,:,nSli];

datatime  = vec(traj_mat["datatime"]);
dt        = traj_mat["tdwelltraj"];
kspha_raw = traj_mat["phs_spha"][:, 1:9] ./ (2π);
StartTime = 50;

nX, nY, nCha = size(csm); nZ = 1
nSample = length(datatime)


##
g = Grid(nX=nX, nY=nY, nZ=nZ, Δx=1.5, Δy=1.5, Δz=1, x=vec(x), y=vec(y), z=vec(z))
BHO = BlochHighOrder("111")
Nblocks = 10;
use_gpu = true;
verbose = true;
######################################################################################
# test with different JumpFacts, delay_applied
######################################################################################
kspha = InterpTrajTime(kspha_raw, dt, StartTime, datatime);
HOOp = HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
data = HOOp * vec(im0);
im1 = HOOp' * data;
data = reshape(data, nSample, nCha);
signalAmpl = sum(abs.(data), dims=1)/ nSample;

iter_max=30
for snr in [Inf, 10, 5];
    JumpFacts = [1,3,6,10]
    τ_true = collect(-5:0.1:5)
    τ_auto = zeros(length(JumpFacts), length(τ_true))
    for i = eachindex(JumpFacts), j = eachindex(τ_true)
        delay_applied = τ_true[j] 
        JumpFact = JumpFacts[i]
        
        if snr == Inf
            data_in = data
        else
            data_in = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));
        end

        τ = FindDelay_v2(g, Complex{T}.(data_in), T.(kspha_raw), T.(datatime), T.(StartTime - delay_applied), T.(dt);
            JumpFact=JumpFact, fieldmap=T.(b0), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks, iter_max=iter_max, Δτ_min=0.0001)
        τ_auto[i, j] = τ
        @info "JumpFact: $(JumpFact), delay_applied: $(delay_applied), delay_auto: $(τ)"
    end
    matwrite("$(path)/20241231/FindDelay_v2_data_itermax$(iter_max)_snr$(snr)_min1e-4.mat", Dict("JumpFacts"=>JumpFacts, "tau_true"=>τ_true, "tau_auto"=>τ_auto))
    err = (Float64.(τ_auto) .- τ_true') .* 1e3; # [ns]
    ys = [err[idx,:] for idx=eachindex(JumpFacts)];
    xs = [τ_true       for idx=eachindex(JumpFacts)];
    labels = ["γ = $(JumpFacts[idx])" for idx=eachindex(JumpFacts)];
    fig = plt_plot(ys; xs=xs, width=7, height=4.5, labels=labels,
        xlabel=L"τ_{true} \ [us]", ylabel=L"τ_{auto} \ - \ τ_{true} \ [ns]",
        fontsize_label=9, fontsize_legend=7, fontsize_ticklabel=7)
    ax = fig.axes[1]
    ax.xaxis.set_ticks(-5:1:5)
    fig.savefig("$(path)/20241231/FindDelay_v2_fig_itermax$(iter_max)_snr$(snr)_min1e-4.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)
end

# width = height = 8
# fig = plt_scatter([Float64.(τ_auto[1,:])]; labels=["γ = 1"],  xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")
# fig = plt_scatter([Float64.(τ_auto[2,:])]; labels=["γ = 3"],  xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")
# fig = plt_scatter([Float64.(τ_auto[3,:])]; labels=["γ = 6"],  xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")
# fig = plt_scatter([Float64.(τ_auto[4,:])]; labels=["γ = 10"], xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")

JumpFacts = [1,3,6,10]

τ_true = matread("$(path)/FindDelay_v2_data_snrInf_min1e-4.mat")["tau_true"]
τ_auto = matread("$(path)/FindDelay_v2_data_snrInf_min1e-4.mat")["tau_auto"]
err = (Float64.(τ_auto) .- τ_true') .* 1e3; # [ns]
ys = [err[idx,:] for idx=eachindex(JumpFacts)];
xs = [τ_true       for idx=eachindex(JumpFacts)];
labels = ["γ = $(JumpFacts[idx])" for idx=eachindex(JumpFacts)];
fig = plt_plot(ys; xs=xs, width=8, height=5, labels=labels,
    xlabel="τ_true [us]", ylabel="τ_auto - τ_true [ns]",
    fontsize_label=9, fontsize_legend=8, fontsize_ticklabel=8)
ax = fig.axes[1]
ax.xaxis.set_ticks(-5:1:5)
# ax.hlines([-5, 5], -5, 5, color="#CCCCCC", linestyle="-", linewidth=0.5)
fig.savefig("$(path)/fig_matmri_tau_snrInf.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)

matplotlib.rc("mathtext", default="regular")
for snr in [5, 10, Inf], Δτ_min in ["5e-3", "1e-3", "1e-4"]
    matfile = "$(path)/20241231/FindDelay_v2_data_itermax30_snr$(snr)_min$(Δτ_min).mat"
    τ_true = matread(matfile)["tau_true"]
    τ_auto = matread(matfile)["tau_auto"]
    err = (Float64.(τ_auto) .- τ_true') .* 1e3; # [ns]
    ys = [err[idx,:] for idx=eachindex(JumpFacts)];
    xs = [τ_true       for idx=eachindex(JumpFacts)];
    labels = ["γ = $(JumpFacts[idx])" for idx=eachindex(JumpFacts)];
    fig = plt_plot(ys; xs=xs, width=7, height=4.5, labels=labels,
        xlabel=L"τ_{true} \ [us]", ylabel=L"τ_{auto} \ - \ τ_{true} \ [ns]",
        fontsize_label=9, fontsize_legend=7, fontsize_ticklabel=7)
    ax = fig.axes[1]
    ax.xaxis.set_ticks(-5:1:5)
    fig.savefig("$(path)/20241231/FindDelay_v2_fig_itermax30_snr$(snr)_min$(Δτ_min).png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)
end