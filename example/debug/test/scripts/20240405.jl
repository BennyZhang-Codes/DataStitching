# simulation
using KomaHighOrder
hoseq = demo_hoseq()
plot_hoseqd(hoseq)

# phantom
obj = brain_phantom2D(brain2D(); ss=3, location=0.8); info(obj);
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π


# scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(ho0=true, ho1=false, ho2=true);
# sim_params["Nblocks"] = 150;
sim_params["gpu"] = true;
sim_params["gpu_device"] = 0;
sim_params["return_type"]="mat";

# simulate
signal = simulate(obj, hoseq, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq, :nominal)
plot_image(reconstruct_2d_image(raw); title="$(sim_params["sim_method"]) Nominal", height=700, width=750)


protocolName = "$(hoseq.SEQ.DEF["Name"])_101_nominal"
raw.params["protocolName"] = protocolName
path = @__DIR__
mrd = ISMRMRDFile(path*"/test/recon/$(protocolName).mrd")
save(mrd, raw)


B1 = 4.92e-6
Trf = 3.2e-3
zmax = 2e-2
fmax = 5e3
z = range(-zmax, zmax, 400)
Gz = fmax / (γ * zmax)

seq = PulseDesigner.RF_sinc(B1, Trf, sys; G=[0;0;Gz], TBP=8)
p2 = plot_seq(seq; max_rf_samples=Inf, slider=false)
sim_params = Dict{String, Any}("Δt_rf" => Trf / length(seq.RF.A[1]))
M = simulate_slice_profile(seq; z, sim_params)





obj = brain_phantom3D(brain3D_03();ss=1, start_end=[1,400]); info(obj);
obj = brain_phantom3D(brain3D_03();ss=3, start_end=[180,220]); info(obj);
obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
obj.y = obj.y*0.7;
obj.x = obj.x*0.7;



# scanner & sim_params
sys = Scanner();
sim_params = KomaMRICore.default_sim_params(); 
sim_params["sim_method"] = BlochHighOrder(ho0=false, ho1=false, ho2=false);
sim_params["Nblocks"] = 20; #sim_params["Nthreads"] = 70;
sim_params["gpu"] = true;
sim_params["return_type"]="mat";

# simulate
signal = simulate(obj, hoseq, sys; sim_params);

raw_nominal = signal_to_raw_data(signal, hoseq, :nominal)
raw_measured = signal_to_raw_data(signal, hoseq, :measured)
plot_image(reconstruct_2d_image(raw_nominal); title="$(sim_params["sim_method"]) Nominal", height=700, width=750)
plot_image(reconstruct_2d_image(raw_measured); title="$(sim_params["sim_method"]) Measured", height=700, width=750)

using Dates
mrd = Sys.iswindows() ? ISMRMRDFile("E:/skope/$(seq.DEF["Name"])_nominal_$(Dates.now()).mrd") : ISMRMRDFile("/home/jyzhang/Desktop/skope/$(seq.DEF["Name"])_nominal_$(Dates.now()).mrd")
save(mrd, raw_nominal)
mrd = Sys.iswindows() ? ISMRMRDFile("E:/skope/$(seq.DEF["Name"])_measured_$(Dates.now()).mrd") : ISMRMRDFile("/home/jyzhang/Desktop/skope/$(seq.DEF["Name"])_measured_$(Dates.now()).mrd")
save(mrd, raw_measured)

raw = RawAcquisitionData(ISMRMRDFile("E:/skope/xw_sp2d-1mm-r1_nominal.mrd"))

acqData = AcquisitionData(raw)
acqData.traj[1].circular = false #Removing circular window
C = maximum(2*abs.(acqData.traj[1].nodes[:]))  #Normalize k-space to -.5 to .5 for NUFFT
acqData.traj[1].nodes = acqData.traj[1].nodes[1:2,:] ./ C
Nx, Ny = raw.params["reconSize"][1:2]
recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)
recParams[:densityWeighting] = true
rec = reconstruction(acqData, recParams)
image3d  = reshape(rec.data, Nx, Ny, :)
image2d = (abs.(image3d) * prod(size(image3d)[1:2]))[:,:,1]
plot_image(image2d; title="$(sim_params["sim_method"])", height=700, width=750)









# =========plot=========
plot_grads_cumtrapz(hoseq;darkmode=true, width=1200, height=250)
plot_seq(hoseq;darkmode=true, width=1200, height=250)

