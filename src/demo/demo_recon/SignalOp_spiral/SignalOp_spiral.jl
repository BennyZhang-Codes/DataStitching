# julia -t 4
# using CUDA
# device!(1) 

using MRIReco, KomaHighOrder, MRICoilSensitivities, PlotlyJS, MAT
dir = "$(@__DIR__)/src/demo/demo_recon/SignalOp_spiral/results"
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
# direct NUFFT
#######################################################################################
# acqData.traj[1] = Trajectory(Float32.(K_nominal_adc'[1:3,:]), acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile, circular=false);

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
@time rec = reconstruction(acqData, recParams);
# p_direct_NUFFTOp = plot_image(abs.(rec.data[:,:]); title="direct NUFFTOp", width=650, height=600)
img_direct_NUFFTOp = abs.(rec.data[:,:]);
# savefig(p_direct_NUFFTOp,       dir*"/direct_NUFFTOp.svg", width=650, height=600,format="svg")
#######################################################################################
# iterative NUFFTOp
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
@time rec = reconstruction(acqData, recParams);
# p_iter_NUFFTOp = plot_image(abs.(rec.data[:,:]); title="iterative NUFFTOp", width=650, height=600)
img_iter_NUFFTOp = abs.(rec.data[:,:]);

#######################################################################################
# iterative SignalOp normalized
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"

Op = SignalOp(shape, acqData.traj[1]; Nblocks=9)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp_normalized = plot_image(abs.(rec.data[:,:]); title="iterative SignalOp normalized", width=650, height=600)
img_iter_SignalOp_normalized = abs.(rec.data[:,:]);

#######################################################################################
# iterative SignalOp 
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"

Op = SignalOp(shape, tr_nominal, 1e-3, 1e-3; Nblocks=9)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
p_iter_SignalOp = plot_image(abs.(rec.data[:,:]); title="iterative SignalOp", width=650, height=600)
img_iter_SignalOp = abs.(rec.data[:,:]);

diff = img_iter_SignalOp - img_iter_SignalOp_normalized
plot_image(diff; title="iter_SignalOp - iter_SignalOp_normalized", width=650, height=600, zmin=minimum(diff), zmax=maximum(diff))

#######################################################################################
# iterative HighOrderOp
#######################################################################################
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
BHO_reco = "000"
Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
recParams[:encodingOps] = reshape([Op], 1,1)
@time rec = reconstruction(acqData, recParams);
# p_iter_HighOrderOp = plot_image(abs.(rec.data[:,:]); title="iterative HighOrderOp", width=650, height=600)
img_iter_HighOrderOp = abs.(rec.data[:,:])

imgs_dict = Dict(
    "img_direct_NUFFTOp"=>img_direct_NUFFTOp, 
    "img_iter_NUFFTOp"=>img_iter_NUFFTOp, 
    "img_iter_SignalOp_normalized"=>img_iter_SignalOp_normalized, 
    "img_iter_SignalOp"=>img_iter_SignalOp, 
    "img_iter_HighOrderOp"=>img_iter_HighOrderOp)

MAT.matwrite(dir*"/SignalOp_Simu_000.mat", imgs_dict; compress=true)

##### plots
imgs = Array{Float32,3}(undef, Nx, Ny, 5)
imgs[:,:, 1] = img_direct_NUFFTOp; imgs[:,:, 2] = img_iter_NUFFTOp; imgs[:,:, 3] = img_iter_SignalOp_normalized; imgs[:,:, 4] = img_iter_SignalOp; imgs[:,:, 5] = img_iter_HighOrderOp
p_imgs = plot_imgs(imgs, ["direct NUFFTOp", "NUFFTOp", "SignalOp NormKspace", "SignalOp", "HighOrderOp"]; title="", width=1100, height=250)
savefig(p_imgs, dir*"/SignalOp_compare_nominalrecon.svg", width=1100, height=250,format="svg")


imgs_errors = Array{Float32,3}(undef, Nx, Ny, 3)
imgs_errors[:,:, 1] = img_iter_SignalOp - img_iter_NUFFTOp; 
imgs_errors[:,:, 2] = img_iter_HighOrderOp - img_iter_NUFFTOp; 
imgs_errors[:,:, 3] = img_iter_HighOrderOp - img_iter_SignalOp;
p_imgs_errors = plot_imgs(imgs_errors, ["SignalOp - NUFFTOp", "HighOrderOp - NUFFTOp", "HighOrderOp - SignalOp"]; title="", width=780, height=250)
savefig(p_imgs_errors, dir*"/SignalOp_compare_nominalrecon_errors.svg", width=780, height=250,format="svg")

p = plot_imgs(img_iter_SignalOp_normalized - img_iter_NUFFTOp , ["SignalOp NormKspace - NUFFTOp",]    ; title="", width=350, height=250)
savefig(p, dir*"/SignalOp_SignalOpNormKspace-NUFFTOp.svg", width=350, height=250,format="svg")

plot_imgs(img_iter_SignalOp - img_iter_NUFFTOp , ["SignalOp - NUFFTOp",]    ; title="", width=500, height=500)
plot_imgs(img_iter_HighOrderOp - img_iter_NUFFTOp , ["HighOrderOp - NUFFTOp",] ; title="", width=500, height=500)
plot_imgs(img_iter_HighOrderOp - img_iter_SignalOp, ["HighOrderOp - SignalOp",]; title="", width=500, height=500)


Δx = Δy = 1e-3 # m

x, y = 1:Nx, 1:Ny
x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1
plot_bgcolor = "rgb(22,26,29)"
grid_color = "rgb(40,52,66)"
p1 = plot(scatter3d(x=x, y=y, z=z,       marker=attr(size=0.7), mode="markers"), Layout(
    paper_bgcolor="rgba(0,0,0,0)", 
    scene=attr(xaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
                yaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
                zaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color)),
    font=attr(family="Times New Roman",color="gray"),
    width=500,height=500))
p2 = plot(scatter3d(x=x*Δx, y=y*Δy, z=z, marker=attr(size=0.7), mode="markers"), Layout(
    paper_bgcolor="rgba(0,0,0,0)",
    scene=attr(xaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
                yaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
                zaxis=attr(backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color)),
    font=attr(family="Times New Roman",color="gray"),
    width=500,height=500))
savefig(p1, dir*"/grid_points_index.svg",            width=500,height=500,format="svg")
savefig(p2, dir*"/grid_points_physicalposition.svg", width=500,height=500,format="svg")
