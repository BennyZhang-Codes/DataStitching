include("SimType.jl")
export SimType

include("rawdata/generate_raw.jl")
export generate_raw

include("load_file.jl")
export load_raw, load_seq, load_dfc

export demo, demo_hoseq, demo_sim

function demo() ::Nothing
    @info "demos:"
    @info "    1. seq => load_seq(; seq=\"spiral\", r=2)"
    @info "    2. raw (mrd)  => load_raw(BlochHighOrder(\"000\"); simtype=SimType(\"B0wo_T2w_ss3\"))"
    @info "    3. dfc_grad => load_dfc(; seq=\"spiral\", r=2)"
    @info "    4. hoseq => demo_hoseq()"
    @info "    5. sim (mat)  => raw, image = demo_sim()"
end


function demo_hoseq(;skope::Bool=true, skope_method::Symbol=:Stitched) ::HO_Sequence
    seq = load_seq(); # skope sequence
    seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
    hoseq = HO_Sequence(seq); # hoseq
    if skope
        GR_skope = load_dfc(;dfc_method=skope_method);
        hoseq.GR_skope[2:4, :] = hoseq.SEQ.GR;
        hoseq.GR_skope[:,8] = GR_skope;
    end
    return hoseq
end


function demo_sim(;
    hoseq = demo_hoseq(),
    obj = brain_phantom2D(BrainPhantom(); ss=3, location=0.8),
    sim_params=Dict{String,Any}(),
    sim_method::BlochHighOrder=BlochHighOrder("000"))
    plot_hoseqd(hoseq)

    # phantom
    info(obj)
    obj.Δw .= obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π

    # scanner & sim_params
    sys = Scanner();
    sim_params = KomaMRICore.default_sim_params(sim_params)
    sim_params["sim_method"] = sim_method;
    sim_params["gpu"] = true;
    sim_params["return_type"]="mat";

    # simulate
    signal = simulate(obj, hoseq, sys; sim_params);
    raw = signal_to_raw_data(signal, hoseq, :nominal)
    image = reconstruct_2d_image(raw)
    return raw, image 
end