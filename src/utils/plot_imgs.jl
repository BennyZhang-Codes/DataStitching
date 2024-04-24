
function plot_imgs(
    imgs, 
    subplot_titles; 
    title="HighOrderOp", 
    width=1300, 
    height=200, 
    annotations=[],
    margin_top=40,
    margin_bottom=0,
    margin_left=0,
    margin_right=100,
    )
    fig = make_subplots(rows=1,cols=length(subplot_titles));
    [add_trace!(fig, heatmap(z=imgs[:,:,idx], transpose=false,coloraxis="coloraxis"); row=1, col=idx) for idx in eachindex(subplot_titles)]
    l = Layout(;
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        title=attr(text=title, y=1, x=0, xanchor= "left", yanchor= "top", yref="container"),
        margin=attr(t=margin_top,l=margin_left,r=margin_right,b=margin_bottom),
        coloraxis=attr(colorscale="Greys",xanchor="left",colorbar=attr(x=1,thickness=20)),
        font=attr(family="Times New Roman",size=16,color="gray"),
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