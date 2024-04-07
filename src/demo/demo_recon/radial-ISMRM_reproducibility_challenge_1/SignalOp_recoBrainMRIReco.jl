using HDF5, DelimitedFiles, BenchmarkTools, RegularizedLeastSquares
using MRIReco, MRISampling, MRICoilSensitivities, ProgressMeter
using KomaMRI, PlotlyJS
trials = parse(Int,get(ENV,"NUM_TRIALS","3"))
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10000
BenchmarkTools.DEFAULT_PARAMETERS.samples = trials

path = @__DIR__
path = path*"/test/recon/radial-ISMRM_reproducibility_challenge_1"
f_sensitivity  = path*"/data/sensitivitiesMRIReco.h5"
filename = path*"/data/rawdata_brain_radial_96proj_12ch.h5";

data = permutedims(h5read(filename, "rawdata"),[3,2,1,4]);
traj = permutedims(h5read(filename, "trajectory"),[3,2,1]);
N = 300;
Nc = 12;
T = Float32;
toeplitz = parse(Int,get(ENV,"TOEPLITZ","0"));
oversamplingFactor = parse(Float64,get(ENV,"OVERSAMPLING","1.25"));

@info "Threads = $(Threads.nthreads()) Toeplitz=$(toeplitz)  OversamplingFactor=$(oversamplingFactor)"   

#############################################
# load data and form Acquisition data object
##############################################
tr = Trajectory(T.(reshape(traj[1:2,:,:],2,:) ./ N), 96, 512, circular=false);
dat = Array{Array{Complex{T},2},3}(undef,1,1,1);
dat[1,1,1] = 1.e8.*reshape(data,:,12);
acqData = AcquisitionData(tr, dat, encodingSize=(N,N));

################################
# generate coil sensitivity maps
################################

if !isfile(f_sensitivity)
  @info "Espirit"
  acqDataCart = regrid(acqData, (N,N); cgnr_iter=3);
  sensitivity = espirit(acqDataCart, (6,6), 30, eigThresh_1=0.02, eigThresh_2=0.98);
  h5write(f_sensitivity, "/sensitivity", sensitivity);
else
  sensitivity = h5read(f_sensitivity, "/sensitivity");
end

##########################
# CG-SENSE
##########################
@info "reference reco"
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 20
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = toeplitz == 1
params[:oversamplingFactor] = oversamplingFactor
params[:senseMaps] = Complex{T}.(reshape(sensitivity, N, N, 1, Nc))

@info "undersampled reco"
rf = [1,2,3,4]
img_cg = Array{ComplexF32,3}(undef,N,N,4)
for i = eachindex(rf)
  @info "r=$(rf[i])"
  global acqDataSub = convertUndersampledData(sample_kspace(acqData, T(rf[i]), "regular"))
  img_cg[:,:,i] = reconstruction(acqDataSub, params).data
end

##############################
# CG-SENSE  SignalOp
##############################
params = Dict{Symbol, Any}()
params[:reco] = "multiCoil"
params[:reconSize] = encodingSize(acqData)
params[:regularization] = "L2"
params[:λ] = T(1.e-2)
params[:iterations] = 20
params[:relTol] = 0.0
params[:solver] = "cgnr"
params[:toeplitz] = toeplitz == 1
params[:oversamplingFactor] = oversamplingFactor
params[:senseMaps] = Complex{T}.(reshape(sensitivity, N, N, 1, Nc));

numContr, numChan = numContrasts(acqData), numChannels(acqData)
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, params)
senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan)
ft = SignalOp((N, N), tr; Nblocks=3, use_gpu=true)
smaps = senseMaps[:,:,1,:]
S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1)
Op = DiagOp(ft, numChan) ∘ S 
params[:encodingOps] = reshape([Op], 1,1);

rf = [1,2,3,4]
img_cg_SignalOp = Array{ComplexF32,3}(undef,N,N,4)
@time for i = eachindex(rf)
	@info "r=$(rf[i])"
	global acqDataSub = convertUndersampledData(sample_kspace(acqData, T(rf[i]), "regular"))
	ft = SignalOp((N, N), acqDataSub.traj[1]; Nblocks=9, use_gpu=true)
	smaps = senseMaps[:,:,1,:]
	S = SensitivityOp(reshape(ComplexF64.(smaps),:,numChan),1)
	Op = DiagOp(ft, numChan) ∘ S 
	params[:encodingOps] = reshape([Op], 1,1);
	# SENSE reconstruction while monitoring error
	img_cg_SignalOp[:,:,i] = reconstruction(acqDataSub, params).data
end

##########################
#  NUFFT
########################## 
params = Dict{Symbol,Any}()
params[:reconSize] = (N, N)
params[:densityWeighting] = true

rf = [1,2,3,4]
img_nufft = Array{ComplexF32,4}(undef,N,N,4,numChannels(acqData))
for i = eachindex(rf)
  @info "r=$(rf[i])"
  global acqDataSub = convertUndersampledData(sample_kspace(acqData, T(rf[i]), "regular"))
  img_nufft[:,:,i,:] = reconstruction(acqDataSub, params).data[:,:,1,1,:,1]
end
img_nufft = sqrt.(sum(img_nufft.^2;dims=4)[:,:,:])

###################################
# plots
###################################
p_nufft = plot_image_radial(img_nufft)
p_cg = plot_image_radial(img_cg)
p_cg_SignalOp = plot_image_radial(img_cg_SignalOp)
p_error = plot_image_radial(img_cg .- img_cg_SignalOp)
savefig(p_cg,          path*"/reco/radial_cg.svg",          width=850, height=200,format="svg")
savefig(p_cg_SignalOp, path*"/reco/radial_cg_SignalOp.svg", width=850, height=200,format="svg")
savefig(p_nufft,       path*"/reco/radial_nufft.svg",       width=850, height=200,format="svg")
savefig(p_error,       path*"/reco/radial_error_cg-cg_SignalOp.svg",       width=850, height=200,format="svg")

function plot_image_radial(img; title="CG-SENSE", width=nothing, height=nothing)
	i1 = heatmap(z=abs.(rotr90(img[:,:,1])), transpose=false,coloraxis="coloraxis")
	i2 = heatmap(z=abs.(rotr90(img[:,:,2])), transpose=false,coloraxis="coloraxis")
	i3 = heatmap(z=abs.(rotr90(img[:,:,3])), transpose=false,coloraxis="coloraxis")
	i4 = heatmap(z=abs.(rotr90(img[:,:,4])), transpose=false,coloraxis="coloraxis")
	fig = make_subplots(rows=1,cols=4);
	add_trace!(fig, i1; row=1, col=1);
	add_trace!(fig, i2; row=1, col=2);
	add_trace!(fig, i3; row=1, col=3);
	add_trace!(fig, i4; row=1, col=4);
	l = Layout(;
		paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
		# title=attr(text=title, y=0.92, x=0.5, xanchor= "center", yanchor= "middle"),
		margin=attr(t=20,l=0,r=0,b=0),
		coloraxis=attr(colorscale="Greys",xanchor="left",colorbar=attr(x=1,thickness=20)),
		font=attr(family="Times New Roman",size=16,color="gray"),
		xaxis=attr(domain=(0.01,0.24),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		xaxis2=attr(domain=(0.26,0.49),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		xaxis3=attr(domain=(0.51,0.74),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		xaxis4=attr(domain=(0.76,0.99),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		yaxis=attr(scaleanchor="x", anchor="x", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		yaxis2=attr(scaleanchor="x2", anchor="x2", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		yaxis3=attr(scaleanchor="x3", anchor="x3", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		yaxis4=attr(scaleanchor="x4", anchor="x4", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
		annotations=[
			attr(text="R = 1",yanchor="bottom",xanchor="center",xref="x1 domain",x=0.5,yref="y1 domain",y=1,showarrow=false,font=attr(size=16)),
			attr(text="R = 2",yanchor="bottom",xanchor="center",xref="x2 domain",x=0.5,yref="y2 domain",y=1,showarrow=false,font=attr(size=16)),
			attr(text="R = 3",yanchor="bottom",xanchor="center",xref="x3 domain",x=0.5,yref="y3 domain",y=1,showarrow=false,font=attr(size=16)),
			attr(text="R = 4",yanchor="bottom",xanchor="center",xref="x4 domain",x=0.5,yref="y4 domain",y=1,showarrow=false,font=attr(size=16)),],
	)
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	p = update(fig; layout=l);
	# savefig(p, "C:/Users/82357/Desktop/dpi300.svg", width=1200, height=350,format="svg")
	return p
end

## trajectory plots
tr1 = convertUndersampledData(sample_kspace(acqData, 1, "regular")).traj[1].nodes
tr2 = convertUndersampledData(sample_kspace(acqData, 2, "regular")).traj[1].nodes
tr3 = convertUndersampledData(sample_kspace(acqData, 3, "regular")).traj[1].nodes
tr4 = convertUndersampledData(sample_kspace(acqData, 4, "regular")).traj[1].nodes
s1 = scatter(x=tr1[1,:], y=tr1[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)
s2 = scatter(x=tr2[1,:], y=tr2[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)
s3 = scatter(x=tr3[1,:], y=tr3[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)
s4 = scatter(x=tr4[1,:], y=tr4[2,:],mode="markers", marker=attr(size=1, color="#EF553B"),showlegend=false)
fig = make_subplots(rows=1,cols=4);
add_trace!(fig, s1; row=1, col=1);
add_trace!(fig, s2; row=1, col=2);
add_trace!(fig, s3; row=1, col=3);
add_trace!(fig, s4; row=1, col=4);
l = Layout(;
	width=700, height=220,
	paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
	margin=attr(t=50,l=0,r=0,b=0),
	coloraxis=attr(colorscale="Greys",xanchor="left",colorbar=attr(x=1,thickness=20)),
	font=attr(family="Times New Roman",size=16,color="gray"),
	xaxis=attr(domain=(0.01,0.24),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	xaxis2=attr(domain=(0.26,0.49),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	xaxis3=attr(domain=(0.51,0.74),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	xaxis4=attr(domain=(0.76,0.99),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	yaxis=attr(scaleanchor="x", anchor="x", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	yaxis2=attr(scaleanchor="x2", anchor="x2", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	yaxis3=attr(scaleanchor="x3", anchor="x3", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	yaxis4=attr(scaleanchor="x4", anchor="x4", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	annotations=[
		attr(text="R = 1",yanchor="middle",xanchor="center",xref="x1 domain",x=0.5,yref="y1 domain",y=1,showarrow=false,font=attr(size=16)),
		attr(text="R = 2",yanchor="middle",xanchor="center",xref="x2 domain",x=0.5,yref="y2 domain",y=1,showarrow=false,font=attr(size=16)),
		attr(text="R = 3",yanchor="middle",xanchor="center",xref="x3 domain",x=0.5,yref="y3 domain",y=1,showarrow=false,font=attr(size=16)),
		attr(text="R = 4",yanchor="middle",xanchor="center",xref="x4 domain",x=0.5,yref="y4 domain",y=1,showarrow=false,font=attr(size=16)),],
)
p = update(fig; layout=l)


## sensitivity plots
s = [sensitivity[:,:,1,i] for i=1:12];
s = [s[1] s[2] s[3] s[4]; s[5] s[6] s[7] s[8]; s[9] s[10] s[11] s[12]];
plot_image(abs.(rotr90(s)); title="sensitivity - espirit")

i1 = heatmap(z=abs.(rotr90(s)), transpose=false,coloraxis="coloraxis")
fig = make_subplots(rows=1,cols=1);
add_trace!(fig, i1; row=1, col=1);
l = Layout(;
	paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
	title=attr(text="sensitivity - espirit", y=0.95, x=0.5, xanchor= "center", yanchor= "top"),
	margin=attr(t=50,l=5,r=5,b=5),
	coloraxis=attr(colorscale="Greys",xanchor="left",colorbar=attr(x=1,thickness=20)),
	font=attr(family="Times New Roman",size=16,color="gray"),
	xaxis=attr(showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
	yaxis=attr(scaleanchor="x", anchor="x", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"),
)
p = update(fig; layout=l)
##############
# debug
##############
reconSize, weights, L_inv, sparseTrafo, reg, normalize, encOps, solvername, senseMaps = MRIReco.setupIterativeReco(acqData, params)
encodingOps = encOps
encDims = ndims(trajectory(acqData))
if encDims!=length(reconSize)
error("reco-dimensionality $(length(reconSize)) and encoding-dimensionality $(encDims) do not match")
end

numContr, numChan, numSl, numRep = numContrasts(acqData), numChannels(acqData), numSlices(acqData), numRepetitions(acqData)
encParams = MRIReco.getEncodingOperatorParams(;params...)

# noise decorrelation
senseMapsUnCorr = decorrelateSenseMaps(L_inv, senseMaps, numChan)

# set sparse trafo in reg
reg[1].params[:sparseTrafo] = sparseTrafo

# solve optimization problem
Ireco = zeros(Complex{T}, prod(reconSize), numSl, numContr, numRep)
l = 1
k = 1
if encodingOps != nothing
	E = encodingOps[:,k]
else
	E = encodingOps_parallel(acqData, reconSize, senseMapsUnCorr; slice=k, encParams...)
end
E = encodingOps[:,k]
j = 1
W = WeightingOp(Complex{T}; weights=weights[j], rep=numChan)
kdata = multiCoilData(acqData, j, k, rep=l) .* repeat(weights[j], numChan)
if !isnothing(L_inv)
	kdata = vec(reshape(kdata, :, numChan) * L_inv')
end

EFull = ∘(W, E[j], isWeighting=true)
EFullᴴEFull = normalOperator(EFull)
solver = createLinearSolver(solvername, EFull; AᴴA=EFullᴴEFull, reg=reg, params...)
I = solve(solver, kdata; params...)

if isCircular( trajectory(acqData, j) )
	circularShutter!(reshape(I, reconSize), 1.0)
end
Ireco[:,k,j,l] = I
Ireco_ = reshape(Ireco, MRIReco.volumeSize(reconSize, numSl)..., numContr, 1,numRep)
plot_image(abs.(rotr90(Ireco_[:,:])); title="iter $(params[:iterations])")



import KomaHighOrder: prod_SignalOp, ctprod_SignalOp
function prod_SignalOp(xm::SubArray{T}, x::Vector{Float64}, y::Vector{Float64}, nodes::Matrix{Float64};
Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
    xm = xm[:,1]
    return prod_SignalOp(xm, x, y, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu)
end
function ctprod_SignalOp(xm::SubArray{T}, x::Vector{Float64}, y::Vector{Float64}, nodes::Matrix{Float64}; 
Nblocks::Int64=1, parts::Vector{UnitRange{Int64}}=[1:size(nodes,2)], use_gpu::Bool=false) where T<:Union{Real,Complex}
    xm = xm[:,1]
    return ctprod_SignalOp(xm, x, y, nodes; Nblocks=Nblocks, parts=parts, use_gpu=use_gpu)
end
