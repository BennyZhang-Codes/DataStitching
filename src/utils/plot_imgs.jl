
function plot_imgs(
    imgs, 
    subplot_titles; 
    title="HighOrderOp",
    thememode=:dark, 
    width=nothing, 
    height=nothing, 
    annotations=[],
    margin_top=40,
    margin_bottom=0,
    margin_left=0,
    margin_right=100,
    )
    bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)

    fig = make_subplots(rows=1,cols=length(subplot_titles));
    [add_trace!(fig, heatmap(z=imgs[:,:,idx], transpose=false,coloraxis="coloraxis"); row=1, col=idx) for idx in eachindex(subplot_titles)]
    l = Layout(;
        paper_bgcolor=bgcolor, plot_bgcolor="rgba(0,0,0,0)",
        title=attr(text=title, y=1, x=0, xanchor= "left", yanchor= "top", yref="container"),
        margin=attr(t=margin_top,l=margin_left,r=margin_right,b=margin_bottom),
        coloraxis=attr(colorscale="Greys",xanchor="left",colorbar=attr(x=1,thickness=20)),
        font=attr(family="Times New Roman",size=16,color=text_color),
    )
    if height !== nothing
        l.height = height + margin_top + margin_bottom
    end
    if width !== nothing
        l.width = width + margin_left + margin_right
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
    for idx in eachindex(subplot_titles)
        push!(annotations,attr(text=subplot_titles[idx],yanchor="bottom",xanchor="center",xref="x$idx domain",x=0.5,yref="y$idx domain",y=1,showarrow=false,font=attr(size=16)))
    end
    push!(l.fields, :annotations=>annotations)
    p = PlotlyJS.update(fig; layout=l)
    return p
end

"""
    p = plot_img(image; height, width, zmin, zmax, thememode, title)

Plots an image matrix.

# Arguments
- `image`: (`::Matrix{Number}`) image matrix

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `zmin`: (`::Real`, `=minimum(abs.(image[:]))`) reference value for minimum color
- `zmax`: (`::Real`, `=maximum(abs.(image[:]))`) reference value for maximum color
- `thememode`: (`::Symbol`, `=:dark`) Symbol to indicate thememode style
- `title`: (`::String`, `=""`) plot title

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the image matrix
"""
function plot_img(
      image;
      height = nothing,
      width = nothing,
      zmin = minimum(image[:]),
      zmax = maximum(image[:]),
      thememode = :dark,
      title = "",
      margin_top=40,
      margin_bottom=0,
      margin_left=0,
      margin_right=0,
  )
	#Layout
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	l = Layout(;
    title=attr(text=title, y=1, x=0, xanchor= "left", yanchor= "top", yref="container"),
    # yaxis_title="y",
    # xaxis_title="x",
    xaxis=attr(constrain="domain",showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)",gridcolor="rgba(0,0,0,0)"),
    yaxis=attr(scaleanchor="x", anchor="x", showticklabels=false,showgrid=false,zerolinecolor="rgba(0,0,0,0)",gridcolor="rgba(0,0,0,0)"),
    margin=attr(t=margin_top,l=margin_left,r=margin_right,b=margin_bottom),
	font_color=text_color,
    modebar=attr(orientation="v",bgcolor=bgcolor,color=text_color,activecolor=plot_bgcolor),
    hovermode="closest",
	paper_bgcolor=bgcolor,
	plot_bgcolor="rgba(0,0,0,0)",
    )
    if height !== nothing
        l.height = height
    end
    if width !== nothing
        l.width = width
    end
	#Plot
	p = heatmap(z=image,transpose=false,zmin=zmin,zmax=zmax,colorscale="Greys")
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "autoScale", "resetScale2d", "pan", "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)
	return KomaMRIPlots.plot_koma(p, l; config)
end