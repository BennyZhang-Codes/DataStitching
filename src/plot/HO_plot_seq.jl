"""
    p = HO_plot_seq(hoseq::HO_Sequence; kwargs...)

Plots a high order sequence struct.

# Arguments
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct

# Keywords
- `width`: (`::Integer`, `=nothing`) plot width
- `height`: (`::Integer`, `=nothing`) plot height
- `slider`: (`::Bool`, `=true`) boolean to indicate whether to display a slider
- `show_seq_blocks`: (`::Bool`, `=false`) boolean to indicate whether to display sequence blocks
- `darkmode`: (`::Bool`, `=false`) boolean to indicate whether to display darkmode style
- `range`: (`::Vector{Real}`, `=[]`) time range to be displayed initially
- `title`: (`::String`, `=""`) plot title
- `max_rf_samples`: (`::Integer`, `=100`) maximum number of RF samples

# Returns
- `p`: (`::PlotlyJS.SyncPlot`) plot of the Sequence struct
"""
function HO_plot_seq(
      hoseq::HO_Sequence;
      width=nothing,
      height=nothing,
      slider=false,
      show_seq_blocks=false,
      darkmode=false,
      max_rf_samples=Inf,
      range=[],
      title="",
      xaxis="x",
      yaxis="y",
      showlegend=true
    )
    # Get the samples of the events in the sequence
    samples = HO_get_samples(hoseq; off_val=Inf, max_rf_samples)

    # Define general params and the vector of plots
    idx = ["Gx" "Gy" "Gz"]
    O = size(hoseq.SEQ.RF, 1)
    p = [scatter() for _ in 1:(3 + 2O + 1 + 9)]

    # For GRADs
    p[1] = scatter(x=samples.gx.t*1e3, y=samples.gx.A*1e3, name=idx[1], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
            xaxis=xaxis, yaxis=yaxis, legendgroup="Gx", showlegend=showlegend, marker=attr(color="#636EFA"))
    p[2] = scatter(x=samples.gy.t*1e3, y=samples.gy.A*1e3, name=idx[2], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
            xaxis=xaxis, yaxis=yaxis, legendgroup="Gy", showlegend=showlegend, marker=attr(color="#EF553B"))
    p[3] = scatter(x=samples.gz.t*1e3, y=samples.gz.A*1e3, name=idx[3], hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
            xaxis=xaxis, yaxis=yaxis, legendgroup="Gz", showlegend=showlegend, marker=attr(color="#00CC96"))

    # For RFs
    for j in 1:O
        rf_phase = angle.(samples.rf.A[:,j])
        rf_phase[samples.rf.A[:,j] .== Inf] .= Inf
        p[2j-1+3] = scatter(x=samples.rf.t*1e3, y=abs.(samples.rf.A[:,j])*1e6, name="|B1|", hovertemplate="(%{x:.4f} ms, %{y:.2f} μT)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="|B1|", showlegend=showlegend, marker=attr(color="#AB63FA"))
        p[2j+3] = scatter(x=samples.rf.t*1e3, y=rf_phase, text=ones(size(samples.rf.t)), name="<B1", hovertemplate="(%{x:.4f} ms, ∠B1: %{y:.4f} rad)", visible="legendonly",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="<B1", showlegend=showlegend, marker=attr(color="#FFA15A"))
    end

    # For ADCs
    p[2O+3+1] = scatter(x=samples.adc.t*1e3, y=samples.adc.A*1.0, name="ADC", hovertemplate="(%{x:.4f} ms, %{y:i})",
                xaxis=xaxis, yaxis=yaxis, legendgroup="ADC", showlegend=showlegend, color=marker=attr(color="#19D3F3"))

    # for skope measured gradients
    p[2O+3+1+1] = scatter(x=samples.h0.t*1e3, y=samples.h0.A*1e3, name="h0", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h0", showlegend=showlegend, marker=attr(color="#5f4690"))
    p[2O+3+1+2] = scatter(x=samples.h1.t*1e3, y=samples.h1.A*1e3, name="h1", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h1", showlegend=showlegend, marker=attr(color="#1d6996"))
    p[2O+3+1+3] = scatter(x=samples.h2.t*1e3, y=samples.h2.A*1e3, name="h2", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h2", showlegend=showlegend, marker=attr(color="#38a6a5"))
    p[2O+3+1+4] = scatter(x=samples.h3.t*1e3, y=samples.h3.A*1e3, name="h3", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h3", showlegend=showlegend, marker=attr(color="#0f8554"))
    p[2O+3+1+5] = scatter(x=samples.h4.t*1e3, y=samples.h4.A*1e3, name="h4", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h4", showlegend=showlegend, marker=attr(color="#73af48"))
    p[2O+3+1+6] = scatter(x=samples.h5.t*1e3, y=samples.h5.A*1e3, name="h5", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h5", showlegend=showlegend, marker=attr(color="#edad08"))
    p[2O+3+1+7] = scatter(x=samples.h6.t*1e3, y=samples.h6.A*1e3, name="h6", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h6", showlegend=showlegend, marker=attr(color="#e17c05"))
    p[2O+3+1+8] = scatter(x=samples.h7.t*1e3, y=samples.h7.A*1e3, name="h7", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h7", showlegend=showlegend, marker=attr(color="#cc503e"))
    p[2O+3+1+9] = scatter(x=samples.h8.t*1e3, y=samples.h8.A*1e3, name="h8", hovertemplate="(%{x:.4f} ms, %{y:.2f} mT/m²)",
                    xaxis=xaxis, yaxis=yaxis, legendgroup="h8", showlegend=showlegend, marker=attr(color="#94346e"))
    # Return the plot
    l, config = KomaMRIPlots.generate_seq_time_layout_config(title, width, height, range, slider, show_seq_blocks, darkmode; T0=KomaMRIBase.get_block_start_times(hoseq.SEQ))
    return KomaMRIPlots.plot_koma(p, l; config)
end