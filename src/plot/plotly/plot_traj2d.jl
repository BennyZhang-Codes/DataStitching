

function plot_traj2d(
	traj::Trajectory;
	title="",
	width=nothing,
	height=nothing,
	thememode=:dark
)
	kx, ky = traj.nodes[1,:], traj.nodes[2,:]
	s = PlotlyJS.scattergl(x=kx, y=ky,mode="markers", marker=attr(size=2))

	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	l = Layout(;title=title, hovermode="closest",plot_bgcolor=plot_bgcolor,paper_bgcolor=bgcolor,font_color=text_color,
		modebar=attr(orientation="h", yanchor="bottom", xanchor="right", y=1, x=0, bgcolor=bgcolor, color=text_color, activecolor=plot_bgcolor),
		legend=attr(orientation="h", yanchor="bottom", xanchor="left", y=1, x=0),
		xaxis=attr(title="kx",showticklabels=true,showgrid=true,gridcolor=grid_color,zerolinecolor=grid_color),
		yaxis=attr(title="ky",scaleanchor="x", anchor="x", showticklabels=true,showgrid=true,gridcolor=grid_color,zerolinecolor=grid_color),
		margin=attr(t=50,l=0,r=0,b=0)
	)
	if height !== nothing l.height = height end
	if width !== nothing l.width = width end

	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan",
								"tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)
	return PlotlyJS.plot(s, l; config)
end
