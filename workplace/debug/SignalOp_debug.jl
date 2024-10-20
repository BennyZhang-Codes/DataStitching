using KomaHighOrder
using MRIReco

# , MRICoilSensitivities, PlotlyJS, MAT, ImageQualityIndexes, ImageDistances

raw = demo_raw("000"; folder="woB0_wT2")
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);
acqData.traj[1].circular = false;

hoseq = demo_hoseq();
_, K_nominal_adc, _, K_skope_adc = get_kspace(hoseq; Δt=1);
t_adc = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);
# times = t_adc .- minimum(t_adc)
times = KomaMRIBase.get_adc_sampling_times(hoseq.SEQ);

T2map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:T2, target_fov=(150, 150), target_resolution=(1,1));
plot_image(T2map; title="T2map", width=650, height=600)
T2map = (T2map.<=46*1e-3) .* Inf .+ T2map;
plot_image(exp.(-0.10 ./ T2map))
B0map = brain_phantom2D_reference(BrainPhantom(); ss=3, location=0.8, key=:Δw_fat, target_fov=(150, 150), target_resolution=(0.2,0.2));

tr = Trajectory(K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);

Nx, Ny = raw.params["reconSize"][1:2];
Nblocks=50; use_gpu=true
# shape = (Nx, Ny);
# Δx = Δy = 1e-3
shape = (Nx*5, Ny*5);
Δx = Δy = 0.2e-3

Nx, Ny = shape
x, y = 1:Nx, 1:Ny
x, y = vec(x .+ y'*0), vec(x*0 .+ y') #grid points
x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1
x, y = x * Δx, y * Δy 

nodes = Float64.(kspaceNodes(tr))
k = size(nodes,2) # number of nodes
Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k

n = k÷Nblocks # number of nodes per block
parts = [n for i=1:Nblocks] # number of nodes per block
parts = [1+n*(i-1):n*i for i=1:Nblocks]
if k%Nblocks!= 0
    push!(parts, n*Nblocks+1:k)
end

out = zeros(ComplexF64, size(x, 1));
nodes = Float64.(kspaceNodes(tr));
times = Float64.(readoutTimes(tr));
xm = acqData.kdata[1];
T2map = vec(T2map);
fieldmap = vec(B0map);

if use_gpu
    out = out |> gpu;
    nodes = nodes |> gpu;
    xm = xm |> gpu;
    x = x |> gpu;
    y = y |> gpu;
    times = times |> gpu;
    fieldmap = fieldmap |> gpu
    T2map = T2map |> gpu;
end

for (block, p) = enumerate(parts)
    kx, ky = @view(nodes[1,p]), @view(nodes[2,p])
    ϕ = (x .* kx') .+ (y .* ky') .+ (fieldmap .* times[p]')
    e = exp.(-2*1im*pi*ϕ )#.+ (times[p]./T2map'*1)')
    # e = exp.(-2*1im*pi*ϕ)
    out +=  conj(e) * xm[p] 
end

if use_gpu
    out = out |> cpu;
    T2map = T2map |> cpu;
end
plot_image(abs.(reshape(out, shape)))



