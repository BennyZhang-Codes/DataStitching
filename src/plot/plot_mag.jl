import KomaMRI.KomaMRICore: Mag


function plot_mag(
	ph::HO_Phantom,
    mag::Mag,
	key::Symbol;
	t0=0,
	height=600,
	width=nothing,
	darkmode=false,
	view_2d=false,
	colorbar=true
)
    @assert key in [:xy_mag, :xy_pha, :xy_real, :xy_imag, :z] "key must be :xy_mag, :xy_pha, :xy_real, xy_imag, or :z"
    @assert size(mag.xy, 1) == size(mag.z, 1) == size(ph.x, 1) "Mag and Phantom must have the same number of spins."

	if key == :xy_mag
		spinstate = abs.(mag.xy)
    elseif key == :xy_pha
        spinstate = angle.(mag.xy)
    elseif key == :xy_real
        spinstate = real.(mag.xy)
    elseif key == :xy_imag
        spinstate = imag.(mag.xy)
    else
        spinstate = mag.z
    end

	cmin_key = minimum(spinstate)
	cmax_key = maximum(spinstate)
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
		h = scatter( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
							mode="markers",
							marker=attr(color=spinstate,
										showscale=colorbar,
										colorscale=colormap,
										colorbar=attr(ticksuffix=unit, title=string(key)),
										cmin=cmin_key,
										cmax=cmax_key,
										size=4
										),
							text=round.(spinstate,digits=4),
							hovertemplate="x: %{x:.1f} cm<br>y: %{y:.1f} cm<br><b>$(string(key))</b>: %{text}$unit<extra></extra>")
		else
		h = scatter3d( x=(ph.x .+ ph.ux(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
								y=(ph.y .+ ph.uy(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
								z=(ph.z .+ ph.uz(ph.x,ph.y,ph.z,t0*1e-3))*1e2,
								mode="markers",
								marker=attr(color=spinstate,
											showscale=colorbar,
											colorscale=colormap,
											colorbar=attr(ticksuffix=unit, title=string(key)),
											cmin=cmin_key,
											cmax=cmax_key,
											size=2
											),
								text=round.(spinstate,digits=4),
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