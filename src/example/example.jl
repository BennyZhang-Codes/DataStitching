include("rawdata/generate_raw.jl")
export generate_raw

include("load_file.jl")
export load_raw, load_seq, load_dfc
export load_dfc_mat
export load_hoseq

export demo_sim




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