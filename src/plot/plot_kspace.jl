import KomaMRI.KomaMRIPlots: plot_kspace

"""
    p = plot_kspace(seq::Sequence; width=nothing, height=nothing, darkmode=false)

Plots the k-space of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the k-space of the Sequence struct

# Examples
```julia-repl
julia> seq_file = joinpath(dirname(pathof(KomaMRI)), "../examples/1.sequences/spiral.seq")

julia> seq = read_seq(seq_file)

julia> plot_kspace(seq)
```
"""
function plot_kspace(
      hoseq::HO_Sequence;
      width=nothing,
      height=nothing,
      darkmode=false
  )
    seq = hoseq.SEQ
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	#Calculations of theoretical k-space
	kspace, kspace_adc = get_kspace(seq; Δt=1) #sim_params["Δt"])
	t_adc = KomaMRIBase.get_adc_sampling_times(seq)
	#Colormap
	c_map = [[t, "hsv($(floor(Int,(1-t)*255)), 100, 50)"] for t=range(0,1;length=10)] # range(s,b,N) only works in Julia 1.7.3
	c = "gray"
	c2_idx = []
	counter = 0
	for s in seq
		if is_ADC_on(s)
			N = s.ADC.N[1]
			append!(c2_idx, counter:N+counter-1)
			counter += N
		end
	end
	c2 = interp_map(c_map, c2_idx ./ maximum(c2_idx))
	#Layout
	mink = minimum(kspace_adc,dims=1)
	maxk = maximum(kspace_adc,dims=1)
	dW = maximum(maxk .- mink, dims=2) * .3
	mink .-= dW
	maxk .+= dW
	#Layout
	l = Layout(;
		paper_bgcolor=bgcolor,
		scene=attr(xaxis=attr(title="kx [m⁻¹]",range=[mink[1],maxk[1]],backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
				   yaxis=attr(title="ky [m⁻¹]",range=[mink[2],maxk[2]],backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
				   zaxis=attr(title="kz [m⁻¹]",range=[mink[3],maxk[3]],backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color)),
		modebar=attr(orientation="h",yanchor="bottom",xanchor="right",y=1,x=0,bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
		legend=attr(orientation="h",yanchor="bottom",xanchor="left",y=1,x=0),
		font_color=text_color,
		scene_camera_eye=attr(x=0, y=0, z=1.7),
		scene_camera_up=attr(x=0, y=1., z=0),
		scene_aspectmode="cube",
		margin=attr(t=0,l=0,r=0))
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	#Plot
	p = [scatter() for j=1:3]
	p[1] = scatter3d(x=kspace[:,1],y=kspace[:,2],z=kspace[:,3],mode="lines",
			line=attr(color=c),name="Trajectory",hoverinfo="skip")
	p[2] = scatter3d(x=kspace_adc[:,1],y=kspace_adc[:,2],z=kspace_adc[:,3],text=round.(t_adc*1e3,digits=3),mode="markers",
			line=attr(color=c2),marker=attr(size=2),name="ADC",hovertemplate="kx: %{x:.1f} m⁻¹<br>ky: %{y:.1f} m⁻¹<br>kz: %{z:.1f} m⁻¹<br><b>t_acq</b>: %{text} ms<extra></extra>")
	p[3] = scatter3d(x=[0],y=[0],z=[0],name="k=0",marker=attr(symbol="cross",size=10,color="red"))
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	return plot_koma(p, l; config)
end

