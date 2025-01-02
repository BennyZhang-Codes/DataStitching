using KomaHighOrder
using MRIReco, MRISimulation
using MAT
# using PyPlot
import RegularizedLeastSquares: SolverInfo
using ImageDistances, ImageQualityIndexes
using CUDA
gpu_idx = 3
CUDA.device!(gpu_idx)

datapath = "$(@__DIR__)/workplace/ReconParams/data";    # data directory
path     = "$(@__DIR__)/workplace/ReconParams/out/MultiChannel"; if ispath(path) == false mkpath(path) end     # output directory

#########################################################################################
# 1. Load sequence and dynamic field data
#     a. Load sequence file (*.seq format from pulseq)
#     b. Load dynamic field data (*.mat format)
#     c. Create a HO_Sequence object as defined in KomaHighOrder.jl by combining the 
#        sequence and dynamic field data.
#########################################################################################

seq_file = "$(datapath)/xw_sp2d_7T-1mm-200-r4.seq"   # *.seq file is the pulseq's sequence file
# under-sampling factor: R = 4
# inplane resolution: 1 mm x 1 mm
# FOV: 200 mm x 200 mm
# readout duration: 29 ms
dfc_file = "$(datapath)/xw_sp2d_7T-1mm-200-r4.mat"   # *.mat file contains the dynamic field data from both stitching method and the standard method.

seq = read_seq(seq_file)[3:end-3]; # read_seq function from KomaMRI.jl, return a struct of Sequence
seq.GR[1,:] = -seq.GR[1,:];        # reverse the sign of the gradient (axis x)
plot_seq(seq)                      # plot_seq function from KomaMRI.jl, plot the Sequence


grad = MAT.matread(dfc_file);
Δt = grad["dt"];
nGradSample, nTerm = size(grad["bfieldStitched"])              
dfcStitched = grad["bfieldStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
dfcStandard = grad["bfieldStandard"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
ntStitched  = grad["nSampleAllSegStitched"]
ntStandard  = grad["nSampleAllSegStandard"]
t = Δt * (nGradSample-1);
GR_dfcStitched = reshape([KomaMRIBase.Grad(dfcStitched[idx,:], t, Δt/2, Δt/2, 0) for idx=1:9], :, 1);
GR_dfcStandard = reshape([KomaMRIBase.Grad(dfcStandard[idx,:], t, Δt/2, Δt/2, 0) for idx=1:9], :, 1);


hoseqStitched = HO_Sequence(seq);                       # hoseq, defined in KomaHighOrder.jl
hoseqStandard = HO_Sequence(seq);                       # hoseq, defined in KomaHighOrder.jl
hoseqStitched.GR_dfc[2:4, :] = hoseqStitched.SEQ.GR;    # copy the 1st-order nominal gradient data from the seq object to the hoseq object
hoseqStandard.GR_dfc[2:4, :] = hoseqStandard.SEQ.GR;    # copy the 1st-order nominal gradient data from the seq object to the hoseq object
hoseqStitched.GR_dfc[:,5] = GR_dfcStitched;             # "5" is the index of the readout block in the spiral sequence
hoseqStandard.GR_dfc[:,5] = GR_dfcStandard;             # "5" is the index of the readout block in the spiral sequence
plot_seq(hoseqStitched)
plot_seq(hoseqStandard)


############################################################################################## 
# Setup
############################################################################################## 
T = Float64;
nX = nY = 200; nZ = 1;  # matrix size for recon
Δx = Δy = 1e-3; Δz = 2e-3;

# settings for Simulation
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
location = 0.65
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
# settings for phantom
csm_type  = :real_32cha;      # all values are 1.0 + 0.0im, single channel
csm_nCoil = 32;         
csm_nRow  = 4;
csm_nCol  = 8;

db0_type  = :quadratic;     
db0_max   = :100.;

# 1. phantom
obj = brain_hophantom2D(phantom; ss=ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= B0 ? obj.Δw : obj.Δw * 0;     # set Δw to 0 if B0=false
obj.T2 .= T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 2. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64";
sim_params["gpu_device"]  = gpu_idx;

# 3. simulate
signal = simulate(obj, hoseqStitched, sys; sim_params);          
raw = signal_to_raw_data(signal, hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)

#### 
data = signal[:,:,1];
nSample, nCha = size(data);
# show the effect of noise
raw = signal_to_raw_data(reshape(AddNoise(data, 1.), (nSample, nCha, 1)), hoseqStitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw, nX=nX, nY=nY);
fig_sos = plt_image(rotl90(sqrt.(sum(img_nufft.^2; dims=3))[:,:,1]))
fig_cha = plt_images(mapslices(rotl90, img_nufft,dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)



# Coil-Sensitivity Map
coil = csm_Real_32cha(217, 181, verbose=true);
csm  = get_center_crop(permutedims(coil, (2,1,3)), nX, nY)[:, :, :];

# csm  = get_center_crop(coil, nX, nY)[end:-1:1, :, :];
fig_csm = plt_images(abs.(csm); dim=3, nRow=csm_nRow, nCol=csm_nCol)


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
_,             _, _, kspha_standard = get_kspace(hoseqStandard; Δt=1);
datatime = KomaMRIBase.get_adc_sampling_times(hoseqStitched.SEQ);
# Add a negative sign, as the phase term used in the model is positive.
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
solver = "cgnr"; regularization = "L2"; iter = 10; λ = 0.;
recParams = Dict{Symbol,Any}()
recParams[:reconSize]      = (nX, nY)
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ]              = λ
recParams[:iterations]     = iter
recParams[:solver]         = solver
recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

recParams[:relTol]         = 0.
recParams[:absTol]         = 0.
# recParams[:ρ] = 0.5

Nblocks = 20;
use_gpu = true;
verbose = false;

HOOp = HighOrderOp(gridding, T.(kspha[:, 1:9]'), T.(datatime); sim_method=BHO, tr_nominal=T.(kspha_nominal'), 
                        Nblocks=Nblocks, fieldmap=T.(b0), csm=Complex{T}.(csm), use_gpu=use_gpu, verbose=verbose);
# @time x = recon_HOOp(HOOp, Complex{T}.(data), Complex{T}.(weight), recParams);
@time x = recon_HOOp(HOOp, Complex{T}.(data), recParams);
plt_image(rotl90(abs.(x)); vmaxp=99.9, title="HighOrderOp, stitched: 111 with ΔB₀")

# solverinfo = recParams[:solverInfo];
# solverinfo.nrmse

fig_xlabel         = "Iteration";
width              = 8
height             = 4.5
fontsize_label     = 9
fontsize_legend    = 7
fontsize_ticklabel = 7
color_label        = "#CCCCCC"
ncols_legend       = 2

solver = "cgnr"; regularization = "L2"; iter = 100;

xs = [collect(1:iter)];
λs = [1, 0.5, 0.1, 5.e-2, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-7, 1e-9];
labels = ["λ=1", "λ=0.5", "λ=0.1", "λ=0.05", "λ=0.01", "λ=1e-3", "λ=1e-4", "λ=1e-5", "λ=1e-7", "λ=0"];
noise_levels = [0, 10, 20, 50];
snrs = 100 ./noise_levels;

for idx_snr in eachindex(snrs)
    convM  = Vector{Vector{Float64}}();
    nrmses = Vector{Vector{Float64}}();
    rmses  = Vector{Vector{Float64}}();
    mses   = Vector{Vector{Float64}}();
    ssims  = Vector{Vector{Float64}}();
    snr = snrs[idx_snr];
    data_in = AddNoise(data, snr);
    for λ in λs
        recParams = Dict{Symbol,Any}()
        recParams[:reconSize]      = (nX, nY)
        recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
        recParams[:λ]              = λ
        recParams[:iterations]     = iter
        recParams[:solver]         = solver
        recParams[:solverInfo]     = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)
        recParams[:relTol]         = 0.
        recParams[:absTol]         = 0.

        @info "solver: $solver, regularization: $regularization, snr: $snr, λ: $λ"
        @time x = recon_HOOp(HOOp, Complex{T}.(data_in), recParams);
        # @time x = recon_HOOp(HOOp, Complex{T}.(data_in), Complex{T}.(weight), recParams);
        
        solverinfo = recParams[:solverInfo];
        x_iters = solverinfo.x_iter[2:end]
        m_mse    = Vector{Float64}();
        m_nrmse  = Vector{Float64}();
        m_rmse   = Vector{Float64}();
        m_ssim   = Vector{Float64}();
        for idx in eachindex(x_iters)
            x_iter = abs.(reshape(x_iters[idx], nX, nY));
            push!(m_nrmse, HO_NRMSE(x_ref, x_iter))
            push!(m_rmse , HO_RMSE(x_ref, x_iter))
            push!(m_mse  , HO_MSE(x_ref, x_iter))
            push!(m_ssim , HO_SSIM(x_ref, x_iter))
        end
        push!(convM,  solverinfo.convMeas[2:end]);
        push!(nrmses, m_nrmse);
        push!(mses,  m_mse);
        push!(rmses, m_rmse);
        push!(ssims, m_ssim);
    end
    prefix = "Spiral_R1_111_$(solver)_$(regularization)_snr$(snr)";
    # prefix = "Spiral_R1_111_densityweight_$(solver)_$(regularization)_snr$(snr)";
    matwrite("$(path)/$(prefix).mat", 
        Dict("convM" => convM, "nrmses" => nrmses, "rmses" => rmses, "mses" => mses, "ssims" => ssims))

    fig =  plt_plot(convM; xlabel=fig_xlabel, ylabel="Convergence", xs=xs, labels=labels, width=width, height=height, 
    fontsize_label=fontsize_label, fontsize_legend=fontsize_legend, fontsize_ticklabel=fontsize_ticklabel, color_label=color_label)
    ax = fig.axes[1]; ax.legend(ncols=ncols_legend, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
    fig.savefig("$(path)/$(prefix)_convM.png", dpi=900, transparent=false, bbox_inches="tight", pad_inches=0.0)

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



solver = "admm"; regularization = "TV"; iter = 50;

xs = [collect(1:iter)];
# λs = [1, 0.5, 0.1, 5.e-2, 1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-7, 1e-9];
# labels = ["λ=1", "λ=0.5", "λ=0.1", "λ=0.05", "λ=0.01", "λ=1e-3", "λ=1e-4", "λ=1e-5", "λ=1e-7", "λ=0"];

λs = [100, 50, 10, 5, 1, 0.5, 0.1, 0.05, 0.01];
labels = ["λ=100", "λ=50", "λ=10", "λ=5", "λ=1", "λ=0.5", "λ=0.1", "λ=0.05", "λ=0.01"];

noise_levels = [20, 50]; # 0, 10, 20, 50
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
        recParams[:solverInfo]     = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)
        recParams[:relTol]         = 0.
        recParams[:absTol]         = 0.

        @info "solver: $solver, regularization: $regularization, snr: $snr, λ: $λ"
        @time x = recon_HOOp(HOOp, Complex{T}.(data_in), recParams);
        # @time x = recon_HOOp(HOOp, Complex{T}.(data_in), Complex{T}.(weight), recParams);
        
        solverinfo = recParams[:solverInfo];
        x_iters = solverinfo.x_iter[2:end]
        m_mse    = Vector{Float64}();
        m_nrmse  = Vector{Float64}();
        m_rmse   = Vector{Float64}();
        m_ssim   = Vector{Float64}();
        for idx in eachindex(x_iters)
            x_iter = abs.(reshape(x_iters[idx], nX, nY));
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
    prefix = "Spiral_R1_111_$(solver)_$(regularization)_snr$(snr)";
    # prefix = "Spiral_R1_111_densityweight_$(solver)_$(regularization)_snr$(snr)";

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



