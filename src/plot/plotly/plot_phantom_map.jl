
"""
    p = plot_phantom_map(obj::HO_Phantom, key::Symbol; kwargs...)

Plots a HO_Phantom map for a specific spin parameter given by `key`.

# Arguments
- `obj`: (`::HO_Phantom`) HO_Phantom struct
- `key`: (`::Symbol`, opts: [`:ρ`, `:T1`, `:T2`, `:T2s`, `:x`, `:y`, `:z`]) symbol for
    displaying different parameters of the HO_Phantom spins

# Keywords
- `t0`: (`::Real`, `=0`, `[ms]`) time to see displacement of the HO_Phantom
- `height`: (`::Integer`, `=600`) plot height
- `width`: (`::Integer`, `=nothing`) plot width
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `view_2d`: (`::Bool`, `=false`) boolean to indicate whether to use a 2D scatter plot
- `colorbar`: (`::Bool`, `=true`) boolean to indicate whether to display a colorbar

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the HO_Phantom map for a specific spin parameter

# References
Colormaps from https://github.com/markgriswold/MRFColormaps
Towards Unified Colormaps for Quantitative MRF Data, Mark Griswold, et al. (2018).

# Examples
```julia-repl
julia> obj2D, obj3D = brain_phantom2D(), brain_phantom3D();

julia> plot_phantom_map(obj2D, :ρ)

julia> plot_phantom_map(obj3D, :ρ)
```
"""
function plot_phantom_map(
      ph::HO_Phantom,
      key::Symbol;
      t0=0,
      height=600,
      width=nothing,
      darkmode=false,
      view_2d=false,
      colorbar=true
  )
	path = @__DIR__
	cmin_key = minimum(getproperty(ph,key))
	cmax_key = maximum(getproperty(ph,key))
	if key == :T1 || key == :T2 || key == :T2s
		cmin_key = 0
		factor = 1e3
		unit = " ms"
		if key  == :T1
			cmax_key = 2500/factor
			colors = MAT.matread(path*"/assets/T1cm.mat")["T1colormap"]
			N, _ = size(colors)
			idx = range(0,1;length=N) #range(0,T,N) works in Julia 1.7
			colormap = [[idx[n], "rgb($(floor(Int,colors[n,1]*255)),$(floor(Int,colors[n,2]*255)),$(floor(Int,colors[n,3]*255)))"] for n=1:N]
		elseif key == :T2 || key == :T2s
			if key == :T2
				cmax_key = 250/factor
			end
    		colors = MAT.matread(path*"/assets/T2cm.mat")["T2colormap"]
			N, _ = size(colors)
			idx = range(0,1;length=N) #range(0,T,N) works in Julia 1.7
			colormap = [[idx[n], "rgb($(floor(Int,colors[n,1]*255)),$(floor(Int,colors[n,2]*255)),$(floor(Int,colors[n,3]*255)))"] for n=1:N]
		end
	elseif key == :x || key == :y || key == :z
		factor = 1e2
		unit = " cm"
		colormap="Greys"
	elseif key == :Δw
		factor = 1/(2π)
		unit = " Hz"
		colormap="Greys"
	else
		factor = 1
		cmin_key = 0
		unit=""
		colormap="Greys"
	end
	cmin_key *= factor
	cmax_key *= factor
	x0 = -maximum(abs.([ph.x ph.y ph.z]))*1e2
    xf =  maximum(abs.([ph.x ph.y ph.z]))*1e2
	#Layout
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = Layout(;title=ph.name*": "*string(key),
		xaxis_title="x",
		yaxis_title="y",
		plot_bgcolor=plot_bgcolor,
		paper_bgcolor=bgcolor,
		xaxis_gridcolor=grid_color,
		yaxis_gridcolor=grid_color,
		xaxis_zerolinecolor=grid_color,
		yaxis_zerolinecolor=grid_color,
		font_color=text_color,
		scene=attr(
			xaxis=attr(title="x",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			yaxis=attr(title="y",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			zaxis=attr(title="z",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			aspectmode="manual",
			aspectratio=attr(x=1,y=1,z=1)),
		margin=attr(t=50,l=0,r=0),
		modebar=attr(orientation="h",bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
		xaxis=attr(constrain="domain"),
		yaxis=attr(scaleanchor="x"),
		hovermode="closest")
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	if view_2d
	h = PlotlyJS.scatter( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
						y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
						mode="markers",
						marker=attr(color=getproperty(ph,key)*factor,
									showscale=colorbar,
									colorscale=colormap,
									colorbar=attr(ticksuffix=unit, title=string(key)),
									cmin=cmin_key,
									cmax=cmax_key,
									size=4
									),
						text=round.(getproperty(ph,key)*factor,digits=4),
						hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
	else
	h = PlotlyJS.scatter3d( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							z=(ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							mode="markers",
							marker=attr(color=getproperty(ph,key)*factor,
										showscale=colorbar,
										colorscale=colormap,
										colorbar=attr(ticksuffix=unit, title=string(key)),
										cmin=cmin_key,
										cmax=cmax_key,
										size=2
										),
							text=round.(getproperty(ph,key)*factor,digits=4),
							hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br>z: %{z:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
	end
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
	)
	return plot_koma(h, l; config)
end


function plot_phantom_map_csm(
	ph::HO_Phantom,
	key::Symbol;
	coil_idx=1,
	t0=0,
	height=600,
	width=nothing,
	darkmode=false,
	view_2d=false,
	colorbar=true
)
	_, nCoil = size(ph.csm)
	@assert coil_idx <= nCoil "Invalid coil index. Must be between 1 and $nCoil"
	if key == :real
		csm = real.(ph.csm[:,coil_idx])
	elseif key == :imag
		csm = imag.(ph.csm[:,coil_idx])
	elseif key == :mag
		csm = abs.(ph.csm[:,coil_idx])
	elseif key == :pha
		csm = angle.(ph.csm[:,coil_idx])
	end

	cmin_key = minimum(csm)
	cmax_key = maximum(csm)
	unit=""
	colormap="Greys"

	x0 = -maximum(abs.([ph.x ph.y ph.z]))*1e2
    xf =  maximum(abs.([ph.x ph.y ph.z]))*1e2
	#Layout
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = theme_chooser(darkmode)
	l = Layout(;title=ph.name*": "*string(key),
		xaxis_title="x",
		yaxis_title="y",
		plot_bgcolor=plot_bgcolor,
		paper_bgcolor=bgcolor,
		xaxis_gridcolor=grid_color,
		yaxis_gridcolor=grid_color,
		xaxis_zerolinecolor=grid_color,
		yaxis_zerolinecolor=grid_color,
		font_color=text_color,
		scene=attr(
			xaxis=attr(title="x",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			yaxis=attr(title="y",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			zaxis=attr(title="z",range=[x0,xf],ticksuffix=" cm",backgroundcolor=plot_bgcolor,gridcolor=grid_color,zerolinecolor=grid_color),
			aspectmode="manual",
			aspectratio=attr(x=1,y=1,z=1)),
		margin=attr(t=50,l=0,r=0),
		modebar=attr(orientation="h",bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
		xaxis=attr(constrain="domain"),
		yaxis=attr(scaleanchor="x"),
		hovermode="closest")
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	if view_2d
		h = PlotlyJS.scatter( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							mode="markers",
							marker=attr(color=csm,
										showscale=colorbar,
										colorscale=colormap,
										colorbar=attr(ticksuffix=unit, title=string(key)),
										cmin=cmin_key,
										cmax=cmax_key,
										size=4
										),
							text=round.(csm,digits=4),
							hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
		else
		h = PlotlyJS.scatter3d( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
								y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
								z=(ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
								mode="markers",
								marker=attr(color=csm,
											showscale=colorbar,
											colorscale=colormap,
											colorbar=attr(ticksuffix=unit, title=string(key)),
											cmin=cmin_key,
											cmax=cmax_key,
											size=2
											),
								text=round.(csm,digits=4),
								hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br>z: %{z:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
		end
		config = PlotConfig(
			displaylogo=false,
			toImageButtonOptions=attr(
				format="svg", # one of png, svg, jpeg, webp
			).fields,
			modeBarButtonsToRemove=["zoom", "pan", "tableRotation", "resetCameraLastSave3d", "orbitRotation", "resetCameraDefault3d"]
		)
		return plot_koma(h, l; config)
end