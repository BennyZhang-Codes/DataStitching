using KomaHighOrder
using MRIReco, MRISimulation
using MAT
using PyPlot
import RegularizedLeastSquares: SolverInfo
using ImageDistances, ImageQualityIndexes
############################################################################################## 
# Setup
############################################################################################## 
simtype = SimType(B0=true, T2=false, ss=5)                       # turn on B0, turn off T2, set phantom subsampling to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
maxOffresonance = 200.                                           # set maximum off-resonance frequency in Hz for quadratic B0 map
Nx = Ny = 150;

dir = "Figures/out"; if ispath(dir) == false mkpath(dir) end     # output directory

# 1. sequence
hoseq_stitched = demo_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = demo_hoseq(dfc_method=:Standard)[4:end]   # :Standard

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=0.8, B0type=:quadratic, maxOffresonance=maxOffresonance)
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseq_stitched, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(img_nufft); title="Sim: $(BHO.name), Δw: [-$maxOffresonance,$maxOffresonance] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(maxOffresonance)_reconNUFFT.svg", width=550,height=500,format="svg")

# ΔB₀ map
B0map = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, target_fov=(150, 150), target_resolution=(1,1),
                                   B0type=:quadratic,key=:Δw, maxOffresonance=maxOffresonance); 
fig_b0map = plt_image(rotl90(B0map), title="B0map [-$maxOffresonance, $maxOffresonance] Hz")
x_ref = brain_phantom2D_reference(phantom; ss=simtype.ss, location=0.8, key=:ρ, target_fov=(150, 150), target_resolution=(1,1));

acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_dfc_adc_stitched = get_kspace(hoseq_stitched; Δt=1);
_, _, _, K_dfc_adc_standard = get_kspace(hoseq_standard; Δt=1);

times = KomaMRIBase.get_adc_sampling_times(hoseq_stitched.SEQ);

tr_nominal          = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_stitched     = Trajectory(K_dfc_adc_stitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_standard     = Trajectory(K_dfc_adc_standard'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

# 1. ideally, with BlochHighOrder("000"), we use nominal trajectory and a null B0map
Op1 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("000"); Nblocks=9, fieldmap=Matrix(B0map).*0, grid=1);
# 4. include both ΔB₀ and all order terms of stitched DFC, with BlochHighOrder("111")
Op4 = HighOrderOp((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("111"); Nblocks=9, fieldmap=Matrix(B0map), grid=1);


# for params search

["admm", "cgnr", "fista", "optista", "pogm", "splitBregman"]
["L2", "L1", "L21", "TV"] # "Positive", "Proj"
[1, 0.5, 0.1, 5.e-2, 1.e-2, 5.e-3, 1.e-3, 1.e-4, 1.e-5, 1.e-7, 1.e-9]

solver = "cgnr"
regularization = "L2"

for solver in ["admm"]
    regularizations = solver == "cgnr" ? ["L2"] : ["L1"]
    for regularization in regularizations
        noise_levels = [0, 2, 4, 6, 8, 10];
        snrs = 100 ./noise_levels;
        convM  = Vector{Vector{Float64}}();
        nrmses = Vector{Vector{Float64}}();
        rmses  = Vector{Vector{Float64}}();
        mses   = Vector{Vector{Float64}}();
        ssims  = Vector{Vector{Float64}}();
        for idx_snr in eachindex(snrs)
            snr = snrs[idx_snr];
            acqData_noise = addNoise(acqData, snr);
            recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
            recParams[:reconSize] = (Nx, Ny)  # 150, 150
            recParams[:densityWeighting] = true
            recParams[:reco] = "standard"
            recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
            recParams[:λ] = 1.e-3
            recParams[:iterations] = 100
            recParams[:solver] = solver
            recParams[:solverInfo] = SolverInfo(vec(ComplexF64.(x_ref)), store_solutions=true);

            recParams[:encodingOps] = reshape([Op4], 1,1);
            @time rec = abs.(reconstruction(acqData_noise, recParams).data[:,:]);
            # plt_image(rotl90(rec); title="HighOrderOp, stitched: 111 with ΔB₀")

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
            for idx in eachindex(snrs)
                ax.plot(curves[i][idx], linewidth=linewidth, label="$(noise_levels[idx]) %")
            end
            ax.set_facecolor(color_facecoler)
            ax.set_xlabel("Iteration", color=color_label, fontsize=fontsize_label)
            ax.set_ylabel(labels[i], color=color_label, fontsize=fontsize_label)
            ax.tick_params(axis="both", color=color_label, labelcolor=color_label, labelsize=fontsize_ticklabel)
            for spine in ax.spines  # "left", "right", "bottom", "top"
                ax.spines[spine].set_color(color_label)
            end
            ax.legend(ncols=1, fontsize=fontsize_legend, labelcolor=color_label, frameon=false, handlelength=1, handletextpad=0.5, labelspacing=0.1, columnspacing=1)
        end
        fig.tight_layout(pad=0.3)
        fig.savefig("$(dir)/addNoise_solver_$(solver)_reg_$(regularization).png", dpi=300, bbox_inches="tight")
    end
end
