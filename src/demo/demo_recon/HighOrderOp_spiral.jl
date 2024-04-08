# julia -t 4
# using CUDA
# device!(1) 

using MRIReco, KomaHighOrder, MRICoilSensitivities, PlotlyJS
dir = @__DIR__
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
# HighOrderOp 
# [1] Simu: 111, Reco: 000
# [2] Simu: 111, Reco: 100
# [3] Simu: 111, Reco: 010
# [4] Simu: 111, Reco: 001
# [5] Simu: 111, Reco: 011
# [6] Simu: 111, Reco: 101
# [7] Simu: 111, Reco: 110
# [8] Simu: 111, Reco: 111
#######################################################################################
BHO_simu = "111";raw = demo_raw(BHO_simu);
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);acqData.traj[1].circular = false;shape = (Nx, Ny);

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
BHO_recos = ["000", "100", "010", "001", "011", "101", "110", "111"]
imgs = Array{ComplexF32,3}(undef, Nx, Ny, length(BHO_recos))
for idx in eachindex(BHO_recos)
    @info "Simu: $(BHO_simu), Reco: $(BHO_recos[idx])"
    BHO_reco = BHO_recos[idx]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    imgs[:,:, idx] = rec.data[:,:]
end

imgs_111 = abs.(imgs)
imgs_111_error = Array{Float32,3}(undef, size(imgs_111))
for idx in eachindex(BHO_recos)
    imgs_111_error[:,:, idx] = imgs_111[:,:, idx] - imgs_111[:,:, end]
end

subplot_titles = ["Reco: $t" for t in BHO_recos]
title="HighOrderOp, Simu: $(BHO_simu)"

p_111       = plot_imgs(imgs_111, subplot_titles; title=title, width=1300, height=200)
p_111_error = plot_imgs(imgs_111_error, subplot_titles; title=title*", error map", width=1300, height=200)

savefig(p_111,       dir*"/HighOrderOp_Simu_111.svg", width=1300, height=200,format="svg")
savefig(p_111_error, dir*"/HighOrderOp_Simu_111_errormap.svg", width=1300, height=200,format="svg")

#######################################################################################
# HighOrderOp 
# [1] Simu: 000, Reco: 000
# [2] Simu: 000, Reco: 100
# [3] Simu: 000, Reco: 010
# [4] Simu: 000, Reco: 001
# [5] Simu: 000, Reco: 011
# [6] Simu: 000, Reco: 101
# [7] Simu: 000, Reco: 110
# [8] Simu: 000, Reco: 111
#######################################################################################
BHO_simu = "000";raw = demo_raw(BHO_simu);
Nx, Ny = raw.params["reconSize"][1:2];
acqData = AcquisitionData(raw);acqData.traj[1].circular = false;shape = (Nx, Ny);

recParams = Dict{Symbol,Any}()
recParams[:reconSize] = (Nx, Ny)  # 150, 150
recParams[:densityWeighting] = true
recParams[:reco] = "standard"
recParams[:regularization] = "L2"
recParams[:λ] = 1.e-2
recParams[:iterations] = 20
recParams[:solver] = "cgnr"
BHO_recos = ["000", "100", "010", "001", "011", "101", "110", "111"]
imgs = Array{ComplexF32,3}(undef, Nx, Ny, length(BHO_recos))
for idx in eachindex(BHO_recos)
    @info "Simu: $(BHO_simu), Reco: $(BHO_recos[idx])"
    BHO_reco = BHO_recos[idx]
    Op = HighOrderOp(shape, tr_nominal, tr_skope, BlochHighOrder(BHO_reco);  Nblocks=9)
    recParams[:encodingOps] = reshape([Op], 1,1)
    @time rec = reconstruction(acqData, recParams);
    imgs[:,:, idx] = rec.data[:,:]
end

imgs_000 = abs.(imgs)
imgs_000_error = Array{Float32,3}(undef, size(imgs_000))
for idx in eachindex(BHO_recos)
    imgs_000_error[:,:, idx] = imgs_000[:,:, idx] - imgs_000[:,:, 1]
end

subplot_titles = ["Reco: $t" for t in BHO_recos]
title="HighOrderOp, Simu: $(BHO_simu)"

p_000       = plot_imgs(imgs_000, subplot_titles; title=title, width=1300, height=200)
p_000_error = plot_imgs(imgs_000_error, subplot_titles; title=title*", error map", width=1300, height=200)

savefig(p_000,       dir*"/HighOrderOp_Simu_000.svg", width=1300, height=200,format="svg")
savefig(p_000_error, dir*"/HighOrderOp_Simu_000_errormap.svg", width=1300, height=200,format="svg")

##########################################################################################################
# plot_imgs
##########################################################################################################
function plot_imgs(imgs, subplot_titles; title="HighOrderOp", width=1200, height=330)
    fig = make_subplots(rows=1,cols=length(subplot_titles));
    [add_trace!(fig, heatmap(z=imgs[:,:,idx], transpose=false,coloraxis="coloraxis"); row=1, col=idx) for idx in eachindex(subplot_titles)]
    l = Layout(;
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        title=attr(text=title, y=1, x=0.5, xanchor= "center", yanchor= "top", yref="container"),
        margin=attr(t=40,l=0,r=100,b=0),
        coloraxis=attr(colorscale="Greys",xanchor="left",colorbar=attr(x=1,thickness=20)),
        font=attr(family="Times New Roman",size=16,color="gray"),
    )
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
    for idx in eachindex(subplot_titles)
        Δdomain = 1/length(subplot_titles)
        spacing = Δdomain/100
        if idx == 1
            push!(l.fields, Symbol("xaxis")=>attr(domain=(Δdomain*(idx-1)+spacing, Δdomain*idx-spacing),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"))
            push!(l.fields, Symbol("yaxis")=>attr(scaleanchor="x", anchor="x", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"))
        else 
            push!(l.fields, Symbol("xaxis$idx")=>attr(domain=(Δdomain*(idx-1)+spacing, Δdomain*idx-spacing),showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"))
            push!(l.fields, Symbol("yaxis$idx")=>attr(scaleanchor="x$idx", anchor="x$idx", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)"))
        end
    end
    annotations=[attr(text=subplot_titles[idx],yanchor="bottom",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=1,showarrow=false,font=attr(size=16)) for idx in eachindex(subplot_titles)]
    push!(l.fields, :annotations=>annotations)
    p = update(fig; layout=l)
    return p
end






