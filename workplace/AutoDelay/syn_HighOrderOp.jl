using KomaHighOrder
using MAT
using PyPlot
using MRIReco
using CUDA

CUDA.device!(3)


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
Nblocks = 2;
use_gpu = true;
verbose = true;

solver = "cgnr"; reg = "L2"; iter = 20; λ = 0.
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver

# fig = plt_image(rotr90(abs.(im0)), vmaxp=99.9)
# fig.savefig("$(path)/im0/fig_im0.png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0)

# for delay_applied in [-5.0, -2.5, -0.5, 0, 0.5, 2.5, 5.0];
delay_applied = 2.0;
# creating data without delay
kspha = InterpTrajTime(kspha_raw, dt, StartTime, datatime);
# plt_kspha(kspha_del, dt)
HOOp = HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
data = HOOp * vec(im0);
im1 = HOOp' * data;
data = reshape(data, nSample, nCha);

# snr = 5;
# signalAmpl = sum(abs.(data), dims=1)/ nSample;
# data = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));

# recon with delay applied trajectories
kspha = InterpTrajTime(kspha_raw, dt, StartTime - delay_applied, datatime);
weight = SampleDensity(kspha[:,2:3]', (nX, nY));
HOOp = HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
# @time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
# fig = plt_image(rotr90(abs.(x1)), vmaxp=99.9)
# fig.savefig("$(path)/im0/fig_img_snrInf_del_$(round(delay_applied, digits=2)).png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0)
# end

τ = FindDelay(g, Complex{T}.(data), T.(kspha_raw), T.(datatime), T.(StartTime-delay_applied), T.(dt);
    JumpFact=3, Δτ_min=0.005, iter_max=20, fieldmap=T.(b0), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)

τ = FindDelay(g, Complex{T}.(data), T.(kspha_raw), T.(datatime), T.(StartTime-delay_applied + τ), T.(dt);
    JumpFact=1, Δτ_min=0.000001, iter_max=20, fieldmap=T.(b0), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks)
######################################################################################
# test with gridding 
######################################################################################
im_grid = zeros(size(im0));
for i in 10:10:(size(im_grid, 1)-1) # 绘制水平网格线
    im_grid[i, :] .= 1
end
for j in 10:10:(size(im_grid, 2)-1) # 绘制垂直网格线
    im_grid[:, j] .= 1
end
fig = plt_image(rotr90(abs.(im_grid)), vmaxp=99.9)
fig.savefig("$(path)/im_grid/fig_im_grid.png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0)
# plt_image(im_grid)
for delay_applied in [-5.0, -2.5, -0.5, 0, 0.5, 2.5, 5.0];
    # delay_applied = 2.;
    # creating data without delay
    kspha = InterpTrajTime(kspha_raw, dt, StartTime, datatime);
    # plt_kspha(kspha_del, dt)
    HOOp = HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    data = HOOp * vec(im_grid);
    im1 = HOOp' * data;
    data = reshape(data, nSample, nCha);
    
    # snr = 5;
    # signalAmpl = sum(abs.(data), dims=1)/ nSample;
    # data = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));
    
    # recon with delay applied trajectories
    kspha = InterpTrajTime(kspha_raw, dt, StartTime - delay_applied, datatime);
    weight = SampleDensity(kspha[:,2:3]', (nX, nY));
    HOOp = HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
    @time x1 = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
    fig = plt_image(rotr90(abs.(x1)), vmaxp=99.9)
    fig.savefig("$(path)/im_grid/fig_im_grid_snrInf_del_$(round(delay_applied, digits=2)).png", dpi=300, transparent=false, bbox_inches="tight", pad_inches=0)
end

######################################################################################
# test with different JumpFacts, delay_applied
######################################################################################
kspha = InterpTrajTime(kspha_raw, dt, StartTime, datatime);
HOOp = HighOrderOp(g, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
data = HOOp * vec(im0);
im1 = HOOp' * data;
data = reshape(data, nSample, nCha);

JumpFacts = [1,3,6,10]
τ_true = collect(-5:0.1:5)
τ_auto = zeros(length(JumpFacts), length(τ_true))
snr = 10;
for i = eachindex(JumpFacts), j = eachindex(τ_true)
    delay_applied = τ_true[j] 
    JumpFact = JumpFacts[i]
    
    signalAmpl = sum(abs.(data), dims=1)/ nSample;
    data_in = data + signalAmpl/snr .* ( randn(size(data))+ 1im*randn(size(data)));
    # data_in = data

    τ = FindDelay1(g, Complex{T}.(data_in), T.(kspha_raw), T.(datatime), T.(StartTime - delay_applied), T.(dt);
        JumpFact=JumpFact, fieldmap=T.(b0), csm=Complex{T}.(csm), sim_method=BHO, Nblocks=Nblocks, iter_max=20)
    τ_auto[i, j] = τ
    @info "JumpFact: $(JumpFact), delay_applied: $(delay_applied), delay_auto: $(τ)"
end



matwrite("$(path)/data_tau_snrInf.mat", Dict("JumpFacts"=>JumpFacts, "tau_true"=>τ_true, "tau_auto"=>τ_auto))



width = height = 8
fig = plt_scatter([Float64.(τ_auto[1,:])]; labels=["γ = 1"],  xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")
fig = plt_scatter([Float64.(τ_auto[2,:])]; labels=["γ = 3"],  xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")
fig = plt_scatter([Float64.(τ_auto[3,:])]; labels=["γ = 6"],  xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")
fig = plt_scatter([Float64.(τ_auto[4,:])]; labels=["γ = 10"], xs=[τ_true], axis_equal=true, width=width, height=height, xlabel="τ_true [us]", ylabel="τ_auto [us]")



τ_true = matread("$(path)/data_tau_snr5.mat")["tau_true"]
τ_auto = matread("$(path)/data_tau_snr5.mat")["tau_auto"]
err = (Float64.(τ_auto) .- τ_true') .* 1e3; # [ns]
ys = [err[idx,:] for idx=eachindex(JumpFacts)];
xs = [τ_true       for idx=eachindex(JumpFacts)];
labels = ["γ = $(JumpFacts[idx])" for idx=eachindex(JumpFacts)];
fig = plt_plot(ys; xs=xs, width=12, height=5, labels=labels,
    xlabel="τ_true [us]", ylabel="τ_auto - τ_true [ns]",
    fontsize_label=10, fontsize_legend=10, fontsize_ticklabel=8)
ax = fig.axes[1]
# ax.hlines([-5, 5], -5, 5, color="#CCCCCC", linestyle="-", linewidth=0.5)
fig.savefig("$(path)/fig_matmri_tau_snrInf.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.05)

#####################
# debug
#####################
delay_applied = -0.3
kspha_del = InterpTrajTime(kspha_raw, dt, StartTime+delay_applied, datatime);
HOOp = HighOrderOp(g, T.(kspha_del[:, 1:9]'), T.(datatime); sim_method=BHO, Nblocks=Nblocks, csm=Complex{T}.(csm), fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
data = HOOp * vec(im0);
im1 = HOOp' * data;
data = reshape(data, nSample, nCha);

gridding    = g                    
data        = Complex{T}.(data) 
kspha       = T.(kspha_raw)    
datatime    = T.(datatime)
StartTime   = T.(StartTime)                        
dt          = T.(dt)              
JumpFact    = 3      
intermode   = AkimaMonotonicInterpolation()
fieldmap    = T.(b0)
csm         = Complex{T}.(csm)
sim_method  = BHO
Nblocks     = Nblocks   
use_gpu     = true   
solver      = "cgnr" 
reg         = "L2"   
iter_max    = 20     
λ           = 0.     
verbose     = false

@assert size(data,2) == size(csm,3) "data and csm must have the same number of coil channels"
@assert size(data,1) == size(datatime, 1) "data and datatime must have the same number of spatial points"

nSample, nCha = size(data);

# 1. Initialize some variables
τ         = 0;
nIter     = 1;
Δτ_prev   = Δτ = Inf;
Δτ_min    = 0.005;        # [us]
τ_perIter = Vector{Float64}();
push!(τ_perIter, τ);
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (gridding.nX, gridding.nY)
recParams[:regularization] = reg
recParams[:λ]              = λ
recParams[:iterations]     = iter_max
recParams[:solver]         = solver


# 2. Compute kspha_dt, which is "dk/dt"
kernel = [1/8 1/4 0 -1/4 -1/8]';
kspha_dt= conv(kspha, kernel)[3:end-2,:]/dt;  # crop both sides to match size of kspha

for τ = -0.320:0.005:-0.280
    JumpFact = 1
    # Interpolate to match datatime
    kspha_τ    = T.(InterpTrajTime(kspha   , dt, τ + StartTime, datatime, intermode=intermode)[1:nSample,1:9]');
    kspha_dt_τ = T.(InterpTrajTime(kspha_dt, dt, τ + StartTime, datatime, intermode=intermode)[1:nSample,1:9]');
    
    weight = SampleDensity(kspha_τ[2:3,:], (gridding.nX, gridding.nY));
    W = WeightingOp(Complex{T}, weights=weight, rep=nCha);

    HOOp    = HighOrderOp(gridding, kspha_τ, datatime; sim_method=sim_method, 
                Nblocks=Nblocks, csm=csm, fieldmap=fieldmap, use_gpu=use_gpu, verbose=verbose);
    HOOp_dt = HighOrderOp(gridding, kspha_τ, datatime; sim_method=sim_method, tr_kspha_dt=kspha_dt_τ, 
                Nblocks=Nblocks, csm=csm, fieldmap=fieldmap, use_gpu=use_gpu, verbose=verbose);
    
    # Update image
    x = recon_HOOp(HOOp, data, weight, recParams)
    if verbose
        plt_image(abs.(x), title="Iteration $nIter", vmaxp=99.9)
    end
    # Update delay
    y1 = vec(data) - HOOp * vec(x); # Y - Aₚxₚ
    y2 = HOOp_dt * vec(x);       # Bₚxₚ
    Δτ = JumpFact * real((W*y2) \ (W*y1));
    # @info (W*y2) \ (W*y1)
    @info "Iteration $nIter, w/  W: Δτ = $(round(Δτ/dt, digits=5)) [us] | τ = $(round(τ/dt, digits=5)) [us] | JumpFact = $JumpFact"
    Δτ = JumpFact * real(y2 \ y1);
    # @info y2 \ y1
    @info "Iteration $nIter, w/o W: Δτ = $(round(Δτ/dt, digits=5)) [us] | τ = $(round(τ/dt, digits=5)) [us] | JumpFact = $JumpFact"
end



###############################################################
# delay estimation for stitching dynamic field monitoring 
###############################################################
τ = -0.3
JumpFact = 1
# Interpolate to match datatime
kspha_τ    = T.(InterpTrajTime(kspha   , dt, τ + StartTime, datatime, intermode=intermode)[1:nSample,1:9]');
kspha_dt_τ = T.(InterpTrajTime(kspha_dt, dt, τ + StartTime, datatime, intermode=intermode)[1:nSample,1:9]');

weight = SampleDensity(kspha_τ[2:3,:], (gridding.nX, gridding.nY));
W = WeightingOp(Complex{T}, weights=weight, rep=nCha);

HOOp    = HighOrderOp(gridding, kspha_τ, datatime; sim_method=sim_method, 
            Nblocks=Nblocks, csm=csm, fieldmap=fieldmap, use_gpu=use_gpu, verbose=verbose);
HOOp_dt = HighOrderOp(gridding, kspha_τ, datatime; sim_method=sim_method, tr_kspha_dt=kspha_dt_τ, 
            Nblocks=Nblocks, csm=csm, fieldmap=fieldmap, use_gpu=use_gpu, verbose=verbose);

# Update image
x = recon_HOOp(HOOp, data, weight, recParams)
if verbose
    plt_image(abs.(x), title="Iteration $nIter", vmaxp=99.9)
end
# Update delay
y1 = vec(data) - HOOp * vec(x); # Y - Aₚxₚ
y2 = HOOp_dt * vec(x);       # Bₚxₚ
Δτ = JumpFact * real((W*y2) \ (W*y1));
# @info (W*y2) \ (W*y1)
@info "Iteration $nIter, w/  W: Δτ = $(round(Δτ/dt, digits=5)) [us] | τ = $(round(τ/dt, digits=5)) [us] | JumpFact = $JumpFact"
Δτ = JumpFact * real(y2 \ y1);
# @info y2 \ y1
@info "Iteration $nIter, w/o W: Δτ = $(round(Δτ/dt, digits=5)) [us] | τ = $(round(τ/dt, digits=5)) [us] | JumpFact = $JumpFact"



r = 1:100
real(vec(reshape(y2, nSample, nCha)[r, :]) \ vec(reshape(y1, nSample, nCha)[r, :]))