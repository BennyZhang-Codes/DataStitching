using PlotlyJS
sampling_params=KomaMRIBase.default_sampling_params()
hoseqd = HO_discretize(hoseq; sampling_params)
WebGL=true
scatter_used = WebGL ? scattergl : scatter
Gx  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gx*1e3, name="Gx", mode="markers+lines", marker_symbol=:circle)
Gy  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gy*1e3, name="Gy", mode="markers+lines", marker_symbol=:circle)
Gz  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gz*1e3, name="Gz", mode="markers+lines", marker_symbol=:circle)
B1  = scatter_used(x=hoseqd.seqd.t*1e3, y=abs.(hoseqd.seqd.B1*1e6), name="|B1|", mode="markers+lines", marker_symbol=:circle)
ADC = scatter_used(x=hoseqd.seqd.t[hoseqd.seqd.ADC]*1e3, y=zeros(sum(hoseqd.seqd.ADC)), name="ADC", mode="markers", marker_symbol=:x)
h0  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h0*1e3, name="h0", mode="markers+lines", marker_symbol=:circle)
h1  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h1*1e3, name="h1", mode="markers+lines", marker_symbol=:circle)
h2  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h2*1e3, name="h2", mode="markers+lines", marker_symbol=:circle)
h3  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h3*1e3, name="h3", mode="markers+lines", marker_symbol=:circle)
h4  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h4*1e3, name="h4", mode="markers+lines", marker_symbol=:circle)
h5  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h5*1e3, name="h5", mode="markers+lines", marker_symbol=:circle)
h6  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h6*1e3, name="h6", mode="markers+lines", marker_symbol=:circle)
h7  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h7*1e3, name="h7", mode="markers+lines", marker_symbol=:circle)
h8  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.h8*1e3, name="h8", mode="markers+lines", marker_symbol=:circle)
# KomaMRIPlots.plot_koma([Gx,Gy,Gz,B1,ADC,h0,h1,h2,h3,h4,h5,h6,h7,h8])

l = Layout(;
    colorway=[
            "#636efa", "#EF553B", "#00cc96",  # Gx, Gy, Gz
            "#AB63FA",                        # |B1|
            # "#FFA15A",                      # <B1
            "#19d3f3",                        # ADC                      #"#FF6692", "#B6E880", "#FF97FF", "#FECB52"
            "#5f4690", "#1d6996", "#38a6a5", "#0f8554", "#73af48", "#edad08", "#e17c05", "#cc503e", "#94346e",  # prism, for skope measured gradients
            #"#6f4070", "#666666"
        ],
    xaxis=attr(ticksuffix=" ms"),
    margin=attr(t=0,l=0,r=0,b=0),
)

config = PlotConfig(
    displaylogo=false,
    toImageButtonOptions=attr(
        format="png", # one of png, svg, jpeg, webp
        width=10*300, 
        height=5*300,
        scale=1,
    ).fields,
    modeBarButtonsToRemove=["zoom", "select2d", "lasso2d", "autoScale", "resetScale2d", "pan",
                            "tableRotation", "resetCameraLastSave", "zoomIn", "zoomOut"]
)

KomaMRIPlots.plot_koma([Gx,Gy,Gz]; config)
HO_plot_hoseqd(hoseq;width=1200,height=200)
open("C:/Users/82357/Desktop/plot_hoseq.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
    end