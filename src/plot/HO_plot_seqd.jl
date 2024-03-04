function HO_plot_seqd(
        hoseq::HO_Sequence; 
        sampling_params=KomaMRIBase.default_sampling_params(),
        WebGL=true
    )
	hoseqd = HO_discretize(hoseq; sampling_params)
    scatter_used = WebGL ? scattergl : scatter
    Gx  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gx*1e3, name="Gx", mode="markers+lines", marker_symbol=:circle)
    Gy  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gy*1e3, name="Gy", mode="markers+lines", marker_symbol=:circle)
    Gz  = scatter_used(x=hoseqd.seqd.t*1e3, y=hoseqd.seqd.Gz*1e3, name="Gz", mode="markers+lines", marker_symbol=:circle)
    B1  = scatter_used(x=hoseqd.seqd.t*1e3, y=abs.(hoseqd.seqd.B1*1e6), name="|B1|", mode="markers+lines", marker_symbol=:circle)
    ADC = scatter_used(x=hoseqd.seqd.t[hoseqd.seqd.ADC]*1e3, y=zeros(sum(hoseqd.seqd.ADC)), name="ADC", mode="markers", marker_symbol=:x)
	KomaMRIPlots.plot_koma([Gx,Gy,Gz,B1,ADC])
end

function HO_plot_hoseqd(
    hoseq::HO_Sequence; 
    sampling_params=KomaMRIBase.default_sampling_params(),
    WebGL=true
    )
	hoseqd = HO_discretize(hoseq; sampling_params)
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
    KomaMRIPlots.plot_koma([Gx,Gy,Gz,B1,ADC,h0,h1,h2,h3,h4,h5,h6,h7,h8])
end



