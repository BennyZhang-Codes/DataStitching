using KomaHighOrder
using MRIReco, MRISimulation
using MAT
# using PyPlot
import RegularizedLeastSquares: SolverInfo
using ImageDistances, ImageQualityIndexes
using CUDA
gpu_idx = 3
CUDA.device!(gpu_idx)

path = "$(@__DIR__)/workplace/ReconParams/out"; if ispath(path) == false mkpath(path) end     # output directory

############################################################################################## 
# Setup
############################################################################################## 
T = Float64;
nX = nY = 150; nZ = 1;  # matrix size for recon
Δx = Δy = 1e-3; Δz = 2e-3;

# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
location = 0.8
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
# settings for phantom
csm_type  = :fan;      # all values are 1.0 + 0.0im, single channel
csm_nCoil = 1;         
csm_nRow  = nothing;
csm_nCol  = nothing;

db0_type  = :quadratic;     
db0_max   = :100.;


# 1. sequence
hoseqStitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseqStandard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64";
sim_params["gpu_device"]  = gpu_idx;

# 4. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);
data = signal[:,:,1];
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(img_nufft); title="Sim: $(BHO.name)")

#### 
nSample, nCha = size(data);
# show the effect of noise
raw = signal_to_raw_data(reshape(AddNoise(data, 1.), (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(img_nufft))


# ΔB₀ map
b0map = brain_phantom2D_reference(phantom, :Δw, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
fig_b0map = plt_B0map(rotl90(b0map))

x_ref = brain_phantom2D_reference(phantom, :ρ, (T.(nX), T.(nY)), (T.(Δx*1e3), T.(Δy*1e3)); location=location, ss=ss);
fig_ref = plt_image(rotl90(x_ref))


#############################################################################
# HighOrderOp, the extended signal model for high-order terms
#############################################################################
# get the k coefficients for the nominal (x,y,z) and the stitching measurement (up to second order)
_, kspha_nominal, _, kspha_stitched = get_kspace(hoseqStitched; Δt=1);
datatime = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);

kspha         = -kspha_stitched;        # Add a negative sign, as the phase term used in the model is positive.
kspha_nominal = -kspha_nominal;
b0            = -b0map;            

nSample, nCha = size(data);
weight = SampleDensity(kspha'[2:3,:], (nX, nY));

# For HighOrderOp
x, y = 1:nX, 1:nY;
# x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
# x, y = x .- nX/2 .- 1, y .- nY/2 .- 1

x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- nX/2 .- 1, y .- nY/2 .- 1

x, y = x * Δx, y * Δy; 
gridding = Grid(nX=nX, nY=nY, nZ=nZ, Δx=Δx, Δy=Δy, Δz=Δz, x=T.(x), y=T.(y), z=T.(z));

solver = "admm"; regularization = "TV"; iter = 5; λ = 1e-4;
# solver = "cgnr"; regularization = "L2"; iter = 10; λ = 5e-1;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

recParams[:relTol] = 0.
# recParams[:ρ] = 0.5

Nblocks = 20;
use_gpu = true;
verbose = false;

HOOp = HighOrderOp(gridding, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), use_gpu=use_gpu, verbose=verbose);
@time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
# @time x = recon_HOOp(HOOp, Complex{T}.(data), recParams);
plt_image(rotl90(abs.(x)); vmaxp=99.9, title="HighOrderOp, stitched: 111 with ΔB₀")

solverinfo = recParams[:solverInfo];
solverinfo.nrmse

HO_SSIM(x_ref, abs.(x))
HO_SSIM(abs.(x), abs.(x_ref))

x2 = HO_img_scale1(x_ref, abs.(x))
plt_image(rotl90(x2))

x3 = HO_img_scale1(x_ref, x)
plt_image(rotl90(abs.(x3)))

plt_image(rotl90(abs.(x)))
plt_image(rotl90(x_ref))

fig_xlabel = "Iteration";
width              = 8
height             = 4.5
fontsize_label     = 9
fontsize_legend    = 7
fontsize_ticklabel = 7
color_label        = "#CCCCCC"
ncols_legend       = 2

# solver = "cgnr"; regularization = "L2"; iter = 100;

# xs = [collect(1:iter)];
# λs = [1, 0.5, 0.1, 5.e-2, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-7, 1e-9];
# labels = ["λ=1", "λ=0.5", "λ=0.1", "λ=0.05", "λ=0.01", "λ=1e-3", "λ=1e-4", "λ=1e-5", "λ=1e-7", "λ=0"];
# noise_levels = [0, 10, 20, 50];
# snrs = 100 ./noise_levels;

# for idx_snr in eachindex(snrs)
#     convM  = Vector{Vector{Float64}}();
#     nrmses = Vector{Vector{Float64}}();
#     rmses  = Vector{Vector{Float64}}();
#     mses   = Vector{Vector{Float64}}();
#     ssims  = Vector{Vector{Float64}}();
#     snr = snrs[idx_snr];
#     data_in = AddNoise(data, snr);
#     for λ in λs
#         recParams = Dict{Symbol,Any}()
#         recParams[:reconSize]      = (nX, nY)
#         recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
#         recParams[:λ]              = λ
#         recParams[:iterations]     = iter
#         recParams[:solver]         = solver
#         recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

#         @info "solver: $solver, regularization: $regularization, snr: $snr, λ: $λ"
#         @time x = recon_HOOp(HOOp, Complex{T}.(data), recParams);
        
#         solverinfo = recParams[:solverInfo];
#         x_iters = solverinfo.x_iter[2:end]
    #     m_mse    = Vector{Float64}();
    #     m_nrmse  = Vector{Float64}();
    #     m_rmse   = Vector{Float64}();
    #     m_ssim   = Vector{Float64}();
    #     for idx in eachindex(x_iters)
    #         x_iter = abs.(reshape(x_iters[idx], 150, 150));
    #         push!(m_nrmse, HO_NRMSE(x_ref, x_iter))
    #         push!(m_rmse , HO_RMSE(x_ref, x_iter))
    #         push!(m_mse  , HO_MSE(x_ref, x_iter))
    #         push!(m_ssim , HO_SSIM(x_ref, x_iter))
    #     end
    #     push!(convM,  solverinfo.convMeas[2:end]);
    #     push!(nrmses, m_nrmse);
    #     push!(mses,  m_mse);
    #     push!(rmses, m_rmse);
    #     push!(ssims, m_ssim);
    # end
    # prefix = "Spiral_R1_111_$(solver)_$(regularization)_snr$(snr)";
    # matwrite("$(path)/$(prefix).mat", 
    #     Dict("convM" => convM, "nrmses" => nrmses, "rmses" => rmses, "mses" => mses, "ssims" => ssims))

#     fig =  plt_plot(convM; xlabel=fig_xlabel, ylabel="Convergence", xs=xs, labels=labels, width=width, height=height, 
#     fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
#     ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
#     fig.savefig("$(path)/$(prefix)_convM.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

#     fig = plt_plot(mses; xlabel=fig_xlabel, ylabel="MSE", xs=xs, labels=labels, width=width, height=height, 
#     fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
#     ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
#     fig.savefig("$(path)/$(prefix)_mse.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

#     fig = plt_plot(rmses; xlabel=fig_xlabel, ylabel="RMSE", xs=xs, labels=labels, width=width, height=height, 
#     fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
#     ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
#     fig.savefig("$(path)/$(prefix)_rmse.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

#     fig = plt_plot(nrmses; xlabel=fig_xlabel, ylabel="NRMSE", xs=xs, labels=labels, width=width, height=height, 
#     fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
#     ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
#     fig.savefig("$(path)/$(prefix)_nrmse.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

#     fig = plt_plot(ssims; xlabel=fig_xlabel, ylabel="SSIM", xs=xs, labels=labels, width=width, height=height, 
#     fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
#     ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)    
#     fig.savefig("$(path)/$(prefix)_ssim.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)
# end



solver = "admm"; regularization = "TV"; iter = 100;

xs = [collect(1:iter)];
λs = [1, 0.5, 0.1, 5.e-2, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-7, 1e-9];
labels = ["λ=1", "λ=0.5", "λ=0.1", "λ=0.05", "λ=0.01", "λ=1e-3", "λ=1e-4", "λ=1e-5", "λ=1e-7", "λ=0"];
noise_levels = [0]; # 0, 10, 20, 50
snrs = 100 ./noise_levels;

for idx_snr in eachindex(snrs)
    convM1  = Vector{Vector{Float64}}();
    convM2  = Vector{Vector{Float64}}();
    nrmses  = Vector{Vector{Float64}}();
    rmses   = Vector{Vector{Float64}}();
    mses    = Vector{Vector{Float64}}();
    ssims   = Vector{Vector{Float64}}();
    snr = snrs[idx_snr];
    data_in = AddNoise(data, snr);
    for λ in λs
        recParams = Dict{Symbol,Any}()
        recParams[:reconSize]      = (nX, nY)
        recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
        recParams[:λ]              = λ
        recParams[:iterations]     = iter
        recParams[:solver]         = solver
        recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

        @info "solver: $solver, regularization: $regularization, snr: $snr, λ: $λ"
        @time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
        
        solverinfo = recParams[:solverInfo];
        x_iters = solverinfo.x_iter[2:end]
        m_mse    = Vector{Float64}();
        m_nrmse  = Vector{Float64}();
        m_rmse   = Vector{Float64}();
        m_ssim   = Vector{Float64}();
        for idx in eachindex(x_iters)
            x_iter = abs.(reshape(x_iters[idx], 150, 150));
            push!(m_nrmse, HO_NRMSE(x_ref, x_iter))
            push!(m_rmse , HO_RMSE(x_ref, x_iter))
            push!(m_mse  , HO_MSE(x_ref, x_iter))
            push!(m_ssim , HO_SSIM(x_ref, x_iter))
        end
        push!(convM1,  solverinfo.convMeas[3:2:end]);
        push!(convM2,  solverinfo.convMeas[4:2:end]);
        push!(nrmses, m_nrmse);
        push!(mses,  m_mse);
        push!(rmses, m_rmse);
        push!(ssims, m_ssim);
    end
    prefix = "Spiral_R1_111_$(solver)_densityweight_$(regularization)_snr$(snr)";
    matwrite("$(path)/$(prefix).mat", 
        Dict("convM1" => convM1, "convM2" => convM2, "nrmses" => nrmses, "rmses" => rmses, "mses" => mses, "ssims" => ssims))

    fig =  plt_plot(convM1; xlabel=fig_xlabel, ylabel="Convergence 1", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
    fig.savefig("$(path)/$(prefix)_convM1.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

    fig =  plt_plot(convM2; xlabel=fig_xlabel, ylabel="Convergence 2", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
    fig.savefig("$(path)/$(prefix)_convM2.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

    fig = plt_plot(mses; xlabel=fig_xlabel, ylabel="MSE", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
    fig.savefig("$(path)/$(prefix)_mse.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

    fig = plt_plot(rmses; xlabel=fig_xlabel, ylabel="RMSE", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
    fig.savefig("$(path)/$(prefix)_rmse.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

    fig = plt_plot(nrmses; xlabel=fig_xlabel, ylabel="NRMSE", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
    fig.savefig("$(path)/$(prefix)_nrmse.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

    fig = plt_plot(ssims; xlabel=fig_xlabel, ylabel="SSIM", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)    
    fig.savefig("$(path)/$(prefix)_ssim.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)
end



