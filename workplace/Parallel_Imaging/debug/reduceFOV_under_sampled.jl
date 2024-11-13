using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances

outpath = "$(@__DIR__)/workplace/Parallel_Imaging/debug/out/"; if ispath(outpath) == false mkpath(outpath) end     # output directory
############################################################################################## 
# Setup
############################################################################################## 
T = Float64;
simtype  = SimType(B0=false, T2=false, ss=12)                     # turn on B0, turn off T2, set phantom subsampling to 5
BHO      = BlochHighOrder("000", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom  = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
location = 0.8;

# settings for phantom
csm_type  = :real_32cha;      # a simulated birdcage coil-sensitivity
csm_nCoil = 32;              # 8-channel
csm_nRow  = 4;
csm_nCol  = 8;

db0_type  = :quadratic;     
db0_max   = :0.;


# reduce FOV but keep the same resolution
resolution = 0.3;
0.3/(0.2/8)
fov_reduce = 8
fov = 150 / fov_reduce
grad_scale = 1 / 7.5
Nx = Ny = Int(round(500 * grad_scale, digits=0));
Sx = Sy = fov / Nx;


# 1. sequence
hoseq_stitched = load_hoseq(dfc_method=:Stitched, r=30)[4:end]   # :Stitched
hoseq_standard = load_hoseq(dfc_method=:Standard, r=30)[4:end]   # :Standard
hoseq_stitched.SEQ.GR[1:2,5] = hoseq_stitched.SEQ.GR[1:2,5]      #* grad_scale
p_seq = plot_seq(hoseq_stitched)
# PlotlyJS.savefig(p_seq,  "$(dir)/R30_hoseq.svg", width=1000, height=400, format="svg")

# 2. phantom
obj = brain_hophantom2D(phantom; ss=simtype.ss, location=location, 
                        csm_type=csm_type, csm_nCoil=csm_nCoil, csm_nRow=csm_nRow, csm_nCol=csm_nCol, 
                        db0_type=db0_type, db0_max=db0_max); 
obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0;     # γ*1.5 T*(-3.45 ppm)*1e-6 * 2π
obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf;   # cancel T2 relaxiation
obj.x = obj.x ./ fov_reduce;
obj.y = obj.y ./ fov_reduce;

# 3. scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params()
sim_params["sim_method"]  = BHO;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. Run simulation
signal = simulate(obj, hoseq_stitched, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
imgs_nufft = recon_2d(raw, Nx=Nx, Ny=Ny);

fig_cha = plt_images(mapslices(rotl90, abs.(imgs_nufft), dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)
fig_sos = plt_image(rotl90(sqrt.(sum(imgs_nufft.^2; dims=3))[:,:,1]))

# fig.savefig("$(dir)/R30_reducefov$(fov_reduce)_NUFFT.png", pad_inches=0, dpi=300, bbox_inches="tight")

# 5. get Coil-sensitivity maps and reference image
coil = csm_Real_32cha(91, 76);
coil = get_center_crop(coil, Nx, Ny);

sensitivity = Array{ComplexF32,4}(undef, Nx, Ny, 1, csm_nCoil);
for c = 1:csm_nCoil
    sensitivity[:,:,1,c] = transpose(coil[:,:,c])
end
p_smap_mag = plt_images(mapslices(rotl90, abs.(sensitivity[:,:,1,:]), dims=[1,2]); dim=3, nRow= csm_nRow, nCol=csm_nCol)
# plt_images(permutedims(abs.(sensitivity[:,:,1,:]), (3, 1,2)))

x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (Sx,Sy);
                                location=0.8, ss=simtype.ss);
plt_image(rotl90(x_ref ), title="reference ρ")


##########################################################################################
# SENSE reconstruction
##########################################################################################
acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

solver = "admm"; regularization = "TV"; iter = 50; λ = 1e-4;
# solver = "cgnr"; regularization = "L2"; iter = 50; λ = 1e-3;


recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "multiCoil"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver
recParams[:toeplitz] = false
recParams[:oversamplingFactor] = 2
recParams[:senseMaps] = Complex{T}.(reshape(sensitivity, Nx, Ny, 1, csm_nCoil));
# recParams[:solverInfo] = SolverInfo(vec(Complex{T}.(x_ref)), store_solutions=true)

@time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
fig = plt_image(rotl90(rec))
fig.savefig("$(dir)/R30_reducefov$(fov_reduce)_$(solver)_$(regularization)_$(iter)_$(λ).png", pad_inches=0, dpi=300, bbox_inches="tight")
