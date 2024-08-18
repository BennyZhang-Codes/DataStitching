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


function demo_hoseq(;dfc::Bool=true, dfc_method::Symbol=:Stitched, seqname::String="demo", r::Int64=1) ::HO_Sequence
    seq = load_seq(;seqname=seqname, r=r); # dfc sequence
    seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
    hoseq = HO_Sequence(seq); # hoseq
    if dfc
        GR_dfc, ntStitched, ntStandard = load_dfc(;dfc_method=dfc_method, seqname=seqname, r=r);
        hoseq.GR_dfc[2:4, :] = hoseq.SEQ.GR;
        hoseq.GR_dfc[:,8] = GR_dfc;
    end
    return hoseq
end


function demo_sim()
    hoseq = demo_hoseq()
    
    sim_params=Dict{String,Any}()
    sim_method::BlochHighOrder=BlochHighOrder("000")
    # plot_hoseqd(hoseq)

    # HO_Phantom
    obj = brain_hophantom2D(BrainPhantom(); ss=5, location=0.8, nCoil=9)
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
    image = recon_2d(raw)
    return raw, image 
end