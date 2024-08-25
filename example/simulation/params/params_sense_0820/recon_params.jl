using KomaHighOrder
using MRIReco
using MAT
using PyPlot
import RegularizedLeastSquares: SolverInfo
using ImageDistances, ImageQualityIndexes
############################################################################################## 
# Setup
############################################################################################## 
R=6
simtype = SimType(B0=false, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
csmtype= :real_32cha
nCoil   = 32; nrows=4; ncols=8;
BHO = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 0.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;

dir = "Figures/params_sense_0820/"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
seq = load_seq(seqname="demo", r=R)[4:end]
hoseq = HO_Sequence(seq)
plot_seq(hoseq)

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance, csmtype=csmtype, nCoil=nCoil)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"
sim_params["Nblocks"]     = 20
sim_params["gpu_device"]  = 2
# 4. simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal; sim_params=copy(sim_params));
imgs_nufft = recon_2d(raw; Nx=Nx, Ny=Ny);

# fig_nufft = plt_image(rotl90(sqrt.(sum(imgs_nufft.^2; dims=3))[:,:,1]); title="Sim: $(BHO.name), R=$(R), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
fig = plt_image(rotl90(sqrt.(sum(imgs_nufft.^2; dims=3))[:,:,1]); width=12/2.54, height=12/2.54,title="")
fig.tight_layout(pad=0, w_pad=0, h_pad=0)
fig.savefig("$(dir)/R$(R)_NUFFT.png", pad_inches=0, dpi=300, bbox_inches="tight")


# CSM 
shape = (Nx, Ny);
T = Float32;
coil = csmtype == :real_32cha ? csm_Real_32cha(217, 181) : csm_Birdcage(217, 181, nCoil, relative_radius=1.5);
coil = get_center_crop(coil, Nx, Ny);


sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, nCoil);
for c = 1:nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
p_smap_mag = plot_imgs_subplots(  abs.(sensitivity[:,:,1,:]), nrows, ncols; title="$(nCoil) coils: Coil Sensitivity (Simulation)")
# plt_images(permutedims(abs.(sensitivity[:,:,1,:]), (3, 1,2)))


# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
fig_b0map = plt_image(rotl90(B0map), title="B0map [-$maxOffresonance, $maxOffresonance] Hz")
x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));



acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_dfc_adc = get_kspace(hoseq; Δt=1);
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);

tr_nominal = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc     = Trajectory(K_dfc_adc'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);


# 4. include both ΔB₀ and all order terms of stitched DFC, with BlochHighOrder("111")
Op4 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc , BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);



# for params search

["admm", "cgnr", "fista", "optista", "pogm", "splitBregman"]
["L2", "L1", "L21", "TV", "Positive", "Proj"]
[1, 0.5, 0.1, 5.e-2, 1.e-2, 5.e-3, 1.e-3, 1.e-4, 1.e-5, 1.e-7, 1.e-9]

solver = "admm"
regularization = "TV"
λs     = [1,]
λ = 1e-4


for solver in ["admm"]
    regularizations = solver == "cgnr" ? ["L2"] : ["TV",]
    for regularization in regularizations
        # λs     = [1.e-2, 1.e-3, 1.e-4, 1.e-5, 1.e-6, 0]
        λs     = [1.e-3, 1.e-4, 1.e-5]
        convM  = Vector{Vector{Float64}}();
        nrmses = Vector{Vector{Float64}}();
        rmses  = Vector{Vector{Float64}}();
        mses   = Vector{Vector{Float64}}();
        ssims  = Vector{Vector{Float64}}();
        for idx_λ in eachindex(λs)
            λ = λs[idx_λ]
            @info "Solver: $(solver), Regularization: $(regularization), λ: $(λ)"
            recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
            recParams[:reconSize] = (Nx, Ny)  # 150, 150
            recParams[:densityWeighting] = true
            recParams[:reco] = "multiCoil"
            recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
            recParams[:λ] = λ
            recParams[:iterations] = 70
            recParams[:solver] = solver
            recParams[:solverInfo] = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);
            
            

            numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
            # reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, recParams);
            # senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan);
            # smaps = senseMaps[:,:,1,:]

            smaps = sensitivity[:,:,1,:]
            S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1)
            Op = DiagOp(Op4, numChan) ∘ S 
            recParams[:senseMaps] = ComplexF64.(reshape(sensitivity, Nx, Ny, 1, nCoil));
            recParams[:encodingOps] = reshape([Op], 1,1);
            @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
            plt_image(rotl90(rec); title="HighOrderOp, stitched: 111 with ΔB₀")


            solverinfo = recParams[:solverInfo];
            x_iters = solverinfo.x_iter[2:end]
            m_mse    = Vector{Float64}();
            m_nrmse  = Vector{Float64}();
            m_rmse   = Vector{Float64}();
            m_ssim   = Vector{Float64}();
            for idx in eachindex(x_iters)
                x_iter = abs.(reshape(x_iters[idx], 150, 150));
                push!(m_nrmse, HO_NRMSE(x_iter, x_ref))
                push!(m_rmse , HO_RMSE(x_iter, x_ref))
                push!(m_mse  , HO_MSE(x_iter, x_ref))
                push!(m_ssim , HO_SSIM(x_iter, x_ref))
            end
            push!(convM,  solverinfo.convMeas[2:end]);
            push!(nrmses, m_nrmse);
            push!(mses,  m_mse);
            push!(rmses, m_rmse);
            push!(ssims, m_ssim);
        end
        figure_width       = 10
        figure_height      = 2
        linewidth          = 0.8
        fontsize_legend    = 6
        fontsize_label     = 10
        fontsize_ticklabel = 6
        fontsize_title     = 10
        color_facecoler    = "#1F1F1F"
        color_label        = "#CCCCCC"


        fig, axs = plt.subplots(1,3, figsize=(figure_width,figure_height), facecolor=color_facecoler)
        fig.suptitle("Solver: $(solver), Regularization: $(regularization)", x=0.5, y=1, ha="center", va="top", color=color_label, fontsize=fontsize_title)
        curves = [convM, nrmses, ssims]
        labels = ["convergenceMeasure", "NRMSE", "SSIM"]
        for i in eachindex(curves)
            ax = axs[i]
            for idx_λ in eachindex(λs)
                ax.plot(curves[i][idx_λ], linewidth=linewidth, label=string(λs[idx_λ]))
            end
            ax.set_facecolor(color_facecoler)
            ax.set_xlabel("Iteration", color=color_label, fontsize=fontsize_label)
            ax.set_ylabel(labels[i], color=color_label, fontsize=fontsize_label)
            ax.tick_params(axis="both", color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
            for spine in ax.spines  # "left", "right", "bottom", "top"
                ax.spines[spine].set_color(color_label)
            end
            ax.legend(ncols=2, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
        end
        fig.tight_layout(pad=0.3)
        fig.savefig("$(dir)/$(solver)_$(regularization)_R$(R).png", dpi=300, bbox_inches="tight")
    end
end

solvers = ["cgnr", "admm", "admm", "admm"]
regs = ["L2", "L2", "L1", "TV"]
λs = [1.e-4, 1e-4, 1e-3, 1e-4]
iters = [30, 30, 20, 20]


for (solver, reg, λ, iter) in zip(solvers, regs, λs, iters)
    @info "Solver: $(solver), Regularization: $(reg), λ: $(λ), iterations: $(iter)"
    recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
    recParams[:reconSize] = (Nx, Ny)  # 150, 150
    recParams[:densityWeighting] = true
    recParams[:reco] = "multiCoil"
    recParams[:regularization] = reg  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
    recParams[:λ] = λ
    recParams[:iterations] = iter
    recParams[:solver] = solver
    recParams[:solverInfo] = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);


    numContr, numChan = MRIReco.numContrasts(acqData), MRIReco.numChannels(acqData);
    smaps = sensitivity[:,:,1,:];
    S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1);
    Op = DiagOp(Op4, numChan) ∘ S ;
    recParams[:senseMaps] = ComplexF64.(reshape(sensitivity, Nx, Ny, 1, nCoil));
    recParams[:encodingOps] = reshape([Op], 1,1);
    @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
    fig = plt_image(rotl90(rec); width=12/2.54, height=12/2.54,title="")
    fig.tight_layout(pad=0, w_pad=0, h_pad=0)
    fig.savefig("$(dir)/R$(R)_$(solver)_$(reg)_$(λ)_$(iter).png", pad_inches=0, dpi=300, bbox_inches="tight")
end
