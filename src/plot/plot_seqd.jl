import KomaMRI.KomaMRIPlots: plot_seqd

function plot_seqd(
        hoseq::HO_Sequence, hoseqd::HO_DiscreteSequence; 
        WebGL=true,
        width=nothing,
        height=nothing,
        slider=false,
        show_seq_blocks=false,
        thememode=:dark,
        range=[],
        title="",
        xaxis="x",
        yaxis="y",
        showlegend=true
)
    scatter_used = WebGL ? scattergl : scatter
    Gx  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gx*1e3, name="Gx", mode="markers+lines", marker_symbol=:circle, 
            xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
    Gy  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gy*1e3, name="Gy", mode="markers+lines", marker_symbol=:circle, 
            xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
    Gz  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gz*1e3, name="Gz", mode="markers+lines", marker_symbol=:circle, 
            xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
    B1  = scatter_used(x=hoseqd.seqd.t*1e3, y=abs.(hoseqd.seqd.B1*1e6), name="|B1|", mode="markers+lines", marker_symbol=:circle, 
            xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} μT)")
    ADC = scatter_used(x=hoseqd.seqd.t[hoseqd.seqd.ADC]*1e3, y=zeros(sum(hoseqd.seqd.ADC)), name="ADC", mode="markers", marker_symbol=:x, 
            xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:i})")
    # Return the plot
    l, config = plot_hoseqd_generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, thememode; T0=KomaMRIBase.get_block_start_times(hoseq.SEQ))
    KomaMRIPlots.plot_koma([Gx,Gy,Gz,B1,ADC], l; config)
end

function plot_seqd(
        hoseq::HO_Sequence; 
        sampling_params=KomaMRIBase.default_sampling_params(),
        WebGL=true,
        width=nothing,
        height=nothing,
        slider=false,
        show_seq_blocks=false,
        thememode=:dark,
        range=[],
        title="",
        xaxis="x",
        yaxis="y",
        showlegend=true
)
    hoseqd = discretize(hoseq; sampling_params)
    plot_seqd(
        hoseq, hoseqd; 
        WebGL=WebGL, 
        width=width, 
        height=height, 
        slider=slider, 
        show_seq_blocks=show_seq_blocks, 
        thememode=thememode, 
        range=range, 
        title=title, 
        xaxis=xaxis, 
        yaxis=yaxis,
        showlegend=showlegend
    )
end

function plot_hoseqd(
        hoseq::HO_Sequence, hoseqd::HO_DiscreteSequence; 
        WebGL=true,
        width=nothing,
        height=nothing,
        slider=false,
        show_seq_blocks=false,
        thememode=:dark,
        range=[],
        title="",
        xaxis="x",
        yaxis="y",
        showlegend=true
)
        scatter_used = WebGL ? scattergl : scatter
        Gx  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gx*1e3, name="Gx", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
        Gy  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gy*1e3, name="Gy", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
        Gz  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gz*1e3, name="Gz", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
        B1  = scatter_used(x=hoseqd.seqd.t*1e3, y=abs.(hoseqd.seqd.B1*1e6), name="|B1|", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} μT)")
        ADC = scatter_used(x=hoseqd.seqd.t[hoseqd.seqd.ADC]*1e3, y=zeros(sum(hoseqd.seqd.ADC)), name="ADC", mode="markers", marker_symbol=:x, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:i})")
        h0  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h0*1e3, name="h0", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT)")
        h1  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h1*1e3, name="h1", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
        h2  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h2*1e3, name="h2", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
        h3  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h3*1e3, name="h3", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)")
        h4  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h4*1e3, name="h4", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)")
        h5  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h5*1e3, name="h5", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)")
        h6  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h6*1e3, name="h6", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)")
        h7  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h7*1e3, name="h7", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)")
        h8  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h8*1e3, name="h8", mode="markers+lines", marker_symbol=:circle, 
                xaxis=xaxis, yaxis=yaxis, showlegend=showlegend, hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)")
        # Return the plot
        l, config = plot_hoseqd_generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, thememode; T0=KomaMRIBase.get_block_start_times(hoseq.SEQ))
        KomaMRIPlots.plot_koma([Gx,Gy,Gz,B1,ADC,h0,h1,h2,h3,h4,h5,h6,h7,h8], l; config)
end

function plot_hoseqd(
        hoseq::HO_Sequence; 
        sampling_params=KomaMRIBase.default_sampling_params(),
        WebGL=true,
        width=nothing,
        height=nothing,
        slider=false,
        show_seq_blocks=false,
        thememode=:dark,
        range=[],
        title="",
        xaxis="x",
        yaxis="y",
        showlegend=true
)
        hoseqd = discretize(hoseq; sampling_params)
        plot_hoseqd(
            hoseq, hoseqd;
            WebGL=WebGL,
            width=width,
            height=height,
            slider=slider,
            show_seq_blocks=show_seq_blocks,
            thememode=thememode,
            range=range,
            title=title,
            xaxis=xaxis,
            yaxis=yaxis,
            showlegend=showlegend
            )
end

function plot_hoseqd_generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, thememode; T0)
	#LAYOUT
	bgcolor, text_color, plot_bgcolor, grid_color, sep_color = HO_theme_chooser(thememode)
	#Shapes
	shapes = []
    N = length(T0)
	if show_seq_blocks
		aux = [line(
			xref="x", yref="paper",
			x0=T0[i]*1e3, y0=0,
			x1=T0[i]*1e3, y1=1,
			line=attr(color=sep_color, width=2),
			) for i = 1:N]
		append!(shapes, aux)
	end
	l = Layout(;title=title, hovermode="closest",
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
				domain=range[:],
				range=range[:],
				rangeslider=attr(visible=slider),
				rangeselector=attr(
					buttons=[
						attr(count=1,
						label="1m",
						step=10,
						stepmode="backward"),
						attr(step="all")
						]
					),
				),
            colorway=[
                "#636efa", "#EF553B", "#00cc96",  # Gx, Gy, Gz
                "#AB63FA",                        # |B1|
                # "#FFA15A",                      # <B1
                "#19d3f3",                        # ADC                      #"#FF6692", "#B6E880", "#FF97FF", "#FECB52"
                "#5f4690", "#1d6996", "#38a6a5", "#0f8554", "#73af48", "#edad08", "#e17c05", "#cc503e", "#94346e",  # prism, for skope measured gradients
                #"#6f4070", "#666666"
            ],
			shapes = shapes,
			margin=attr(t=0,l=0,r=0,b=0)
		)
	if height !== nothing
		l.height = height
	end
	if width !== nothing
		l.width = width
	end
	#CONFIG
	config = PlotConfig(
		displaylogo=false,
		toImageButtonOptions=attr(
			format="svg", # one of png, svg, jpeg, webp
		).fields,
		modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan",
								"tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
	)

	l, config
end

