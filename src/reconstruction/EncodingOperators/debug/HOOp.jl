using KomaHighOrder
using MRIReco
using MAT
import RegularizedLeastSquares: SolverInfo
using ImageDistances
using ProgressMeter

outpath = "$(@__DIR__)/Figures/Fig5/out"; if ispath(outpath) == false mkpath(outpath) end     # output directory
############################################################################################## 
# Setup
############################################################################################## 
B0 = true     # turn on B0
T2 = false    # turn off T2
ss = 5        # set phantom down-sample factor to 5
BHO = BlochHighOrder("111", true, true)                          # turn on all order terms of dynamic field change, turn on Δw_excitation, Δw_precession
phantom = BrainPhantom(prefix="brain3D724", x=0.2, y=0.2, z=0.2) # decide which phantom file to use
location = 0.8;
Nx = Ny = 150;

# setting the coil sensitivity used in the simulation
csm_type  = :fan;      # a simulated birdcage coil-sensitivity
csm_nCoil = 1;         # 1-channel
csm_nRow  = 1;
csm_nCol  = 1;

db0_type  = :quadratic;     
db0_max   = :200.;            # set the maximum off-resonance frequency in Hz for quadratic B0 map


solver = "admm"; regularization = "TV"; λ = 1.e-4; iter=20;

# 1. sequence
hoseq_stitched = load_hoseq(dfc_method=:Stitched)[4:end]   # :Stitched
hoseq_standard = load_hoseq(dfc_method=:Standard)[4:end]   # :Standard

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
sim_params["gpu"] = true;
sim_params["return_type"] = "mat";
sim_params["precision"]   = "f64"

# 4. simulate
signal = simulate(obj, hoseq_stitched, sys; sim_params);
raw = signal_to_raw_data(signal, hoseq_stitched, :nominal; sim_params=copy(sim_params));
img_nufft = recon_2d(raw);
fig_nufft = plt_image(rotl90(img_nufft); title="Sim: $(BHO.name), Δw: [-$db0_max,$db0_max] Hz")
# savefig(p_image, dir*"/quadraticB0map_$(db0_max)_reconNUFFT.svg", width=550,height=500,format="svg")

# ΔB₀ map (the same as the one used for simulation), we will use this map in reconstruction
B0map = brain_phantom2D_reference(phantom, :Δw, (150., 150.), (1., 1.); location=location, ss=ss, db0_type=db0_type, db0_max=db0_max);
fig_b0map = plt_B0map(rotl90(B0map))

# Proton-density map (reference)
x_ref = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=ss);
fig_ref = plt_image(rotl90(x_ref))

# 5. reconstruction
acqData = AcquisitionData(raw, BHO; sim_params=sim_params);
acqData.traj[1].circular = false;

_, K_nominal_adc, _, K_dfc_adc_stitched = get_kspace(hoseq_stitched; Δt=1);
_, _, _, K_dfc_adc_standard = get_kspace(hoseq_standard; Δt=1);

times = KomaMRIBase.get_adc_sampling_times(hoseq_stitched.SEQ);

tr_nominal          = Trajectory(   K_nominal_adc'[1:3,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_stitched     = Trajectory(K_dfc_adc_stitched'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);
tr_dfc_standard     = Trajectory(K_dfc_adc_standard'[:,:], acqData.traj[1].numProfiles, acqData.traj[1].numSamplingPerProfile; circular=false, times=times);


solver = "cgnr"; regularization = "L2"; λ = 1.e-4; iter=1;
recParams = Dict{Symbol,Any}(); #recParams = merge(defaultRecoParams(), recParams)
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = regularization  # ["L2", "L1", "L21", "TV", "LLR", "Positive", "Proj", "Nuclear"]
recParams[:λ] = λ
recParams[:iterations] = iter
recParams[:solver] = solver

Profile.clear()
@profile begin
    Op = HighOrderOp_i1((Nx, Ny), tr_nominal, tr_dfc_stitched , BlochHighOrder("100"); Nblocks=4, fieldmap=Matrix(B0map), grid=1, use_gpu=true, verbose=true);
    recParams[:encodingOps] = reshape([Op], 1,1);
    @time rec = abs.(reconstruction(acqData, recParams).data[:,:]);
    plt_image(rotl90(rec))
end
Profile.print()


shape        = (Nx, Ny)
tr_nominal   = tr_nominal
tr_measured  = tr_dfc_stitched
sim_method   = BlochHighOrder("000")
fieldmap     = Matrix(B0map)
Nblocks      = 20
use_gpu      = true
verbose      = true
Δx           =1e-3
Δy           =1e-3
grid         =1


nodes_measured = Float64.(kspaceNodes(tr_measured))
nodes_nominal = Float64.(kspaceNodes(tr_nominal))
times = Float64.(readoutTimes(tr_measured))
@assert size(nodes_measured,1) == 9 "nodes for measured must have 9 rows"
@assert size(nodes_nominal,1) == 3 "nodes for nominal must have 3 rows"
@assert size(fieldmap) == shape "fieldmap must have same size as shape"

fieldmap = vec(fieldmap)
nrow = size(nodes_measured,2)
ncol = prod(shape)
k = size(nodes_measured,2) # number of nodes
Nblocks = Nblocks > k ? k : Nblocks # Nblocks must be <= k

n = k÷Nblocks # number of nodes per block
parts = [n for i=1:Nblocks] # number of nodes per block
parts = [1+n*(i-1):n*i for i=1:Nblocks]
if k%Nblocks!= 0
    push!(parts, n*Nblocks+1:k)
end

Nx, Ny = shape
x, y = 1:Nx, 1:Ny
if grid == 1      # x up->down, y left->right
    x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- Nx/2 .- 1, y .- Ny/2 .- 1
elseif grid == 2
    x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
elseif grid == 3  # x left->right, y down->up
    x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
elseif grid == 4  # x left->right, y up->down
    x, y, z = vec(x*0.0 .+ y'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
elseif grid == 5  # x right->left, y up->down
    x, y, z = vec(x*0.0 .+ y[end:-1:1]'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
    x, y = x .- (Nx+1)/2, y .- (Ny+1)/2
end
# print(x)
x, y = x * Δx, y * Δy 

xm = brain_phantom2D_reference(phantom, :ρ, (150., 150.), (1., 1.); location=location, ss=ss);
xm = vec(xm)

xm = Vector(xm)
if verbose
    @info "HighOrderOp prod Nblocks=$Nblocks, use_gpu=$use_gpu"
end
out = zeros(ComplexF64,size(nodes_measured,2))
x0 = ones(Float64, size(x))
if use_gpu
    out = out |> gpu
    nodes_measured = nodes_measured |> gpu
    nodes_nominal = nodes_nominal |> gpu
    xm = xm |> gpu
    x0 = x0 |> gpu
    x = x |> gpu
    y = y |> gpu
    z = z |> gpu
    times = times |> gpu
    fieldmap = fieldmap |> gpu
end
progress_bar = Progress(Nblocks)
for (block, p) = enumerate(parts)
    h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                   @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
    hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
    ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
    ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
    ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
            h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
    ϕB0 = times[p] .* fieldmap'
    ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
    @time e = exp.(-2*1im*pi*ϕ)
    out[p] =  e * xm
    if verbose
        next!(progress_bar, showvalues=[(:Nblocks, block)])
    end
end
if use_gpu
    out = out |> cpu
end



e_ϕ = zeros(ComplexF64, size(nodes_measured,2), size(xm, 1))
progress_bar = Progress(Nblocks)
for (block, p) = enumerate(parts)
    h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                   @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
    hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
    ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
    ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
    ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
            h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
    ϕB0 = times[p] .* fieldmap'
    ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
    @time e[p, :] = exp.(-2*1im*pi*ϕ)
    if verbose
        next!(progress_bar, showvalues=[(:Nblocks, block)])
    end
end
if use_gpu
    out = out |> cpu
end


progress_bar = Progress(Nblocks)
for (block, p) = enumerate(parts)
    @time out[p] = e[p, :] * xm
    if verbose
        next!(progress_bar, showvalues=[(:Nblocks, block)])
    end
end
if use_gpu
    out = out |> cpu
end


e = @view e_ϕ[:, :]

@info "HighOrderOp_li calculaiton of the exponentials of total phase, prod Nblocks=$Nblocks, use_gpu=$use_gpu"

x0 = ones(Float64, ncol)
if use_gpu
    println("using GPU")
    nodes_measured = nodes_measured |> gpu
    nodes_nominal = nodes_nominal |> gpu
    x0 = x0 |> gpu
    x = x |> gpu
    y = y |> gpu
    z = z |> gpu
    times = times |> gpu
    fieldmap = fieldmap |> gpu
end


e_ϕ = zeros(ComplexF32, nrow, ncol)
@time begin
    progress_bar = Progress(Nblocks)
    for (block, p) = enumerate(parts)
        h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
                                    @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
        hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
        ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
        ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
        ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
                h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
        ϕB0 = times[p] .* fieldmap'
        ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
        exp_ϕ = exp.(-2*1im*pi*ϕ)
        if use_gpu
            exp_ϕ = exp_ϕ |> cpu
        end
        e[p, :] =  exp_ϕ
        if verbose
            next!(progress_bar, showvalues=[(:Nblocks, block)])
        end
    end
end

e[p, :] |> cpu



p = 1:4400
h0, h1, h2, h3, h4, h5, h6, h7, h8 = @view(nodes_measured[1,p]), @view(nodes_measured[2,p]), @view(nodes_measured[3,p]), @view(nodes_measured[4,p]),
    @view(nodes_measured[5,p]), @view(nodes_measured[6,p]), @view(nodes_measured[7,p]), @view(nodes_measured[8,p]), @view(nodes_measured[9,p])
hx, hy, hz = @view(nodes_nominal[1,p]), @view(nodes_nominal[2,p]), @view(nodes_nominal[3,p])
ϕ0 = sim_method.ho0 ? h0 .* x0' : 0
ϕ1 = sim_method.ho1 ? (h1 .* x') .+ (h2 .* y') .+ (h3 .* z') : (hx .* x') .+ (hy .* y') .+ (hz .* z')
ϕ2 = sim_method.ho2 ? h4 .* (x .* y)' .+ h5 .* (z .* y)' .+ h6 .* (3z.^2-(x.^2 .+ y.^2 .+ z.^2))' .+
    h7 .* (x .* z)' .+ h8 .* (x.^2 .- y.^2)' : 0
ϕB0 = times[p] .* fieldmap'
ϕ = ϕ0 .+ ϕ1 .+ ϕ2 .+ ϕB0
@time e = exp.(-2*1im*pi*ϕ)
out[p] =  e * xm


@time begin
    for (block, p) = enumerate(parts)
        @time a =  e[p, :]|> cpu
    end
end








csm = zeros(2, 3, 4)

for i in 1:4
    csm[:, :, i] = [1 2 3; 4 5 6] .+ i*10
end

