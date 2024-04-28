import KomaMRI.KomaMRIPlots: plot_kspace, interp_map, theme_chooser, plot_koma



"""
    p = plot_kspace(seq::Sequence; width=nothing, height=nothing, thememode=:dark)

Plots the k-space of a Sequence struct.

# Arguments
- `seq`: (`::Sequence`) Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `thememode`: (`::Bool`, `=false`) boolean to indicate whether to display thememode style

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
      thememode=:dark
  )
    seq = hoseq.SEQ
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	#Calculations of theoretical k-space
	K_nominal, K_nominal_adc, K_skope, K_skope_adc = get_kspace(hoseq; Δt=1) #sim_params["Δt"])
	K_skope = K_skope[:, 2:4]          # H1, H2, H3 => x, y, z
	K_skope_adc = K_skope_adc[:, 2:4]

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
	mink = minimum(K_nominal_adc,dims=1)
	maxk = maximum(K_nominal_adc,dims=1)
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
	p = [scatter() for j=1:5]
	p[1] = scatter3d(x=K_nominal[:,1],y=K_nominal[:,2],z=K_nominal[:,3],mode="lines",
			line=attr(color=c),name="nominal Traj",hoverinfo="skip")
	p[2] = scatter3d(x=K_nominal_adc[:,1],y=K_nominal_adc[:,2],z=K_nominal_adc[:,3],text=round.(t_adc*1e3,digits=3),mode="markers",
			line=attr(color=c2),marker=attr(size=2),name="nominal ADC",hovertemplate="nominal<br>kx: %{x:.1f} m⁻¹<br>ky: %{y:.1f} m⁻¹<br>kz: %{z:.1f} m⁻¹<br><b>t_acq</b>: %{text} ms<extra></extra>")
	p[3] = scatter3d(x=K_skope[:,1],y=K_skope[:,2],z=K_skope[:,3],mode="lines",
			line=attr(color=c),name="measured Traj",hoverinfo="skip")
	p[4] = scatter3d(x=K_skope_adc[:,1],y=K_skope_adc[:,2],z=K_skope_adc[:,3],text=round.(t_adc*1e3,digits=3),mode="markers",
			line=attr(color=c2),marker=attr(size=2),name="measured ADC",hovertemplate="measured<br>kx: %{x:.1f} m⁻¹<br>ky: %{y:.1f} m⁻¹<br>kz: %{z:.1f} m⁻¹<br><b>t_acq</b>: %{text} ms<extra></extra>")
	
	p[5] = scatter3d(x=[0],y=[0],z=[0],name="k=0",marker=attr(symbol="cross",size=10,color="red"))
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	return plot_koma(p, l; config)
end

function plot_kspace(
	hoseq::HO_Sequence,
	key::Symbol;
	width=nothing,
	height=nothing,
	thememode=:dark
)
	@assert key == :x || key == :y || key == :z || key == :all "key must be one of :x, :y, :z, :all"
	seq = hoseq.SEQ
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	K_nominal, K_nominal_adc, K_skope, K_skope_adc = get_kspace(hoseq; Δt=1)
	K_skope = K_skope[:, 2:4]          # H1, H2, H3 => x, y, z
	K_skope_adc = K_skope_adc[:, 2:4]

	t_adc = KomaMRIBase.get_adc_sampling_times(seq)

	t, Δt = KomaMRIBase.get_variable_times(hoseq.SEQ; Δt=1)
	t = t[1:end-1]
	t_seq = t .+ Δt

	#Layout
	mink = minimum(K_nominal_adc,dims=1)
	maxk = maximum(K_nominal_adc,dims=1)
	dW = maximum(maxk .- mink, dims=2) * .3
	mink .-= dW
	maxk .+= dW
	#Layout
	l = Layout(;hovermode="closest",
		xaxis_title="",
		modebar=attr(orientation="h", yanchor="bottom", xanchor="right", y=1, x=0, bgcolor=bgcolor, color=text_color, activecolor=plot_bgcolor),
		legend=attr(orientation="h", yanchor="bottom", xanchor="left", y=1, x=0),
		plot_bgcolor=plot_bgcolor,
		paper_bgcolor=bgcolor,
		xaxis_gridcolor=grid_color,
		yaxis_gridcolor=grid_color,
		xaxis_zerolinecolor=grid_color,
		yaxis_zerolinecolor=grid_color,
		font_color=text_color,
		yaxis_fixedrange = false,
		xaxis=attr(
			ticksuffix=" ms",
			),
		margin=attr(t=0,l=0,r=0,b=0)
	)
	if height !== nothing
	l.height = height
	end
	if width !== nothing
	l.width = width
	end
	#Plot
	px = [scattergl() for j=1:5]
	px[1] = scattergl(x=t_seq*1e3, y=K_nominal[:,1],mode="lines", line=attr(color="#636efa"),name="x",hoverinfo="skip",legendgroup="nominal", legendgrouptitle_text="nominal",hovertemplate="nominal<br>kx: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	px[2] = scattergl(x=t_seq*1e3, y=K_skope[:,1],mode="lines", line=attr(color="#EF553B"),name="x",hoverinfo="skip",legendgroup="measured", legendgrouptitle_text="measured",hovertemplate="measured<br>kx: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	px[3] = scattergl(x=t_seq*1e3, y=K_nominal[:,1]-K_skope[:,1],mode="lines", line=attr(color="#00cc96"),name="x",hoverinfo="skip",legendgroup="diff", legendgrouptitle_text="difference",hovertemplate="diff<br>kx: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	px[4] = scattergl(x=t_adc*1e3, y=K_nominal_adc[:,1],mode="markers",line=attr(color="#19d3f3"),marker=attr(size=5, symbol=:x),name="x",legendgroup="nominal ADC",legendgrouptitle_text="nominal ADC",hovertemplate="nominal ADC<br>kx: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	px[5] = scattergl(x=t_adc*1e3, y=K_skope_adc[:,1],mode="markers",line=attr(color="#FFA15A"),marker=attr(size=5, symbol=:x),name="x",legendgroup="measured ADC",legendgrouptitle_text="measured ADC",hovertemplate="measured ADC<br>kx: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	
	py = [scattergl() for j=1:5]
	py[1] = scattergl(x=t_seq*1e3, y=K_nominal[:,2],mode="lines", line=attr(color="#636efa"),name="y",hoverinfo="skip",legendgroup="nominal",legendgrouptitle_text="nominal",hovertemplate="nominal<br>ky: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	py[2] = scattergl(x=t_seq*1e3, y=K_skope[:,2],mode="lines", line=attr(color="#EF553B"),name="y",hoverinfo="skip",legendgroup="measured",legendgrouptitle_text="measured",hovertemplate="measured<br>ky: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	py[3] = scattergl(x=t_seq*1e3, y=K_nominal[:,2]-K_skope[:,2],mode="lines", line=attr(color="#00cc96"),name="y",hoverinfo="skip",legendgroup="diff",legendgrouptitle_text="difference",hovertemplate="diff<br>ky: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	py[4] = scattergl(x=t_adc*1e3, y=K_nominal_adc[:,2],mode="markers",line=attr(color="#19d3f3"),marker=attr(size=5, symbol=:x),name="y",legendgroup="nominal ADC",legendgrouptitle_text="nominal ADC",hovertemplate="nominal ADC<br>ky: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	py[5] = scattergl(x=t_adc*1e3, y=K_skope_adc[:,2],mode="markers",line=attr(color="#FFA15A"),marker=attr(size=5, symbol=:x),name="y",legendgroup="measured ADC",legendgrouptitle_text="measured ADC",hovertemplate="measured ADC<br>ky: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	
	pz = [scattergl() for j=1:5]
	pz[1] = scattergl(x=t_seq*1e3, y=K_nominal[:,3],mode="lines", line=attr(color="#636efa"),name="z",hoverinfo="skip",legendgroup="nominal",legendgrouptitle_text="nominal",hovertemplate="nominal<br>kz: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	pz[2] = scattergl(x=t_seq*1e3, y=K_skope[:,3],mode="lines", line=attr(color="#EF553B"),name="z",hoverinfo="skip",legendgroup="measured",legendgrouptitle_text="measured",hovertemplate="measured<br>kz: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	pz[3] = scattergl(x=t_seq*1e3, y=K_nominal[:,3]-K_skope[:,3],mode="lines", line=attr(color="#00cc96"),name="z",hoverinfo="skip",legendgroup="diff",legendgrouptitle_text="difference",hovertemplate="diff<br>kz: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	pz[4] = scattergl(x=t_adc*1e3, y=K_nominal_adc[:,3],mode="markers",line=attr(color="#19d3f3"),marker=attr(size=5, symbol=:x),name="z",legendgroup="nominal ADC",legendgrouptitle_text="nominal ADC",hovertemplate="nominal ADC<br>kz: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	pz[5] = scattergl(x=t_adc*1e3, y=K_skope_adc[:,3],mode="markers",line=attr(color="#FFA15A"),marker=attr(size=5, symbol=:x),name="z",legendgroup="measured ADC",legendgrouptitle_text="measured ADC",hovertemplate="measured ADC<br>kz: %{y:.1f} m⁻¹<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	if key == :x
		p = plot_koma(px, l; config)
	elseif key == :y
		p = plot_koma(py, l; config)
	elseif key == :z
		p = plot_koma(pz, l; config)
	elseif key == :all
		p = [plot_koma(px, l; config); plot_koma(py, l; config); plot_koma(pz, l; config)]
	end

	return p
end


function plot_grads_cumtrapz(	
	hoseq::HO_Sequence,
	order::Integer;
	width=nothing,
	height=nothing,
	thememode=:dark)
	@assert 0 <= order <= 8 "order must be between 0 and 8"
	seq = hoseq.SEQ
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	_, _, K_skope, K_skope_adc = get_kspace(hoseq; Δt=1)
	K_skope = K_skope[:, order+1:order+1]          # H1, H2, H3 => x, y, z
	K_skope_adc = K_skope_adc[:, order+1:order+1]

	t_adc = KomaMRIBase.get_adc_sampling_times(seq)

	t, Δt = KomaMRIBase.get_variable_times(hoseq.SEQ; Δt=1)
	t = t[1:end-1]
	t_seq = t .+ Δt

	#Layout
	l = Layout(;hovermode="closest",
		xaxis_title="",
		modebar=attr(orientation="h", yanchor="bottom", xanchor="right", y=1, x=0, bgcolor=bgcolor, color=text_color, activecolor=plot_bgcolor),
		legend=attr(orientation="h", yanchor="bottom", xanchor="left", y=1, x=0),
		plot_bgcolor=plot_bgcolor,
		paper_bgcolor=bgcolor,
		xaxis_gridcolor=grid_color,
		yaxis_gridcolor=grid_color,
		xaxis_zerolinecolor=grid_color,
		yaxis_zerolinecolor=grid_color,
		font_color=text_color,
		yaxis_fixedrange = false,
		xaxis=attr(
			ticksuffix=" ms",
			),
		margin=attr(t=0,l=0,r=0,b=0)
	)
	if height !== nothing
	l.height = height
	end
	if width !== nothing
	l.width = width
	end
	SH = SphericalHarmonics()
	name = SH.dict["h$(order)"].name
	unit = SH.dict["h$(order)"].cumtrapz_unit
	expression = SH.dict["h$(order)"].expression
	l.xaxis_title = "t_seq (ms)"
	l.yaxis_title = "$(name) ($(expression) $(unit))"
	#Plot

	p = [scattergl() for j=1:9]
    p[1] = scattergl(x=t_seq*1e3, y=K_skope[:,1],mode="lines", line=attr(color="#636efa"),name="h$(order)",hoverinfo="skip",hovertemplate="$(name)[$(expression)]: %{y:.3f} $(unit)<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	p[2] = scattergl(x=t_adc*1e3, y=K_skope_adc[:,1],mode="markers",line=attr(color="#19d3f3"),marker=attr(size=5, symbol=:x),name="h$(order) ADC",hovertemplate="$(name)[$(expression)]: %{y:.3f} $(unit)<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	p = plot_koma(p, l; config)
	return p
end


function plot_grads_cumtrapz(	
	hoseq::HO_Sequence;
	width=nothing,
	height=nothing,
	thememode=:dark)
	seq = hoseq.SEQ
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	_, _, K_skope, K_skope_adc = get_kspace(hoseq; Δt=1)

	t_adc = KomaMRIBase.get_adc_sampling_times(seq)
	t, Δt = KomaMRIBase.get_variable_times(hoseq.SEQ; Δt=1)
	t = t[1:end-1]
	t_seq = t .+ Δt
	colors = ["#5f4690" "#1d6996" "#38a6a5" "#0f8554" "#73af48" "#edad08" "#e17c05" "#cc503e" "#94346e"] # prism, for skope measured gradients
	#Layout
	l = Layout(;hovermode="closest",
		xaxis_title="",
		modebar=attr(orientation="h", yanchor="bottom", xanchor="right", y=1, x=0, bgcolor=bgcolor, color=text_color, activecolor=plot_bgcolor),
		legend=attr(orientation="h", yanchor="bottom", xanchor="left", y=1, x=0),
		plot_bgcolor=plot_bgcolor,
		paper_bgcolor=bgcolor,
		xaxis_gridcolor=grid_color,
		yaxis_gridcolor=grid_color,
		xaxis_zerolinecolor=grid_color,
		yaxis_zerolinecolor=grid_color,
		font_color=text_color,
		yaxis_fixedrange = false,
		xaxis=attr(
			ticksuffix=" ms",
			),
		margin=attr(t=0,l=0,r=0,b=0),
		colorway=vec([colors; colors]) # prism, for skope measured gradients
	)
	if height !== nothing
	l.height = height
	end
	if width !== nothing
	l.width = width
	end
	SH = SphericalHarmonics()
	l.xaxis_title = "t_seq (ms)"
	#Plot
	p = [scattergl() for j=1:18]

	for h =0:8
		name = SH.dict["h$(h)"].name
		unit = SH.dict["h$(h)"].cumtrapz_unit
		expression = SH.dict["h$(h)"].expression
		p[1+h*2] = scattergl(x=t_seq*1e3, y=K_skope[:,h+1],mode="lines",name="h$(h)",legendgroup="h$(h)", hoverinfo="skip",hovertemplate="$(name)[$(expression)]: %{y:.3f} $(unit)<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
		p[2+h*2] = scattergl(x=t_adc*1e3, y=K_skope_adc[:,h+1],mode="markers",marker=attr(size=5, symbol=:circle),name="h$(h) ADC",showlegend=false, legendgroup="h$(h)", hovertemplate="h$(h) ADC<br>$(name)[$(expression)]: %{y:.3f} $(unit)<br><b>t_acq</b>: %{x:.3f} ms<extra></extra>")
	end
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	p = plot_koma(p, l; config)
	return p
end