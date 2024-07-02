
"""
generate_raw(SimType(B0=false, T2=false, ss=3))
generate_raw(SimType(B0=false, T2=true,  ss=3))
generate_raw(SimType(B0=false, T2=false, ss=5))
generate_raw(SimType(B0=false, T2=true,  ss=5))
"""

function generate_raw(simtype::SimType)
    folder = simtype.name
    path = "$(@__DIR__)/$(folder)"
    if ispath(path) == false mkdir(path) end

    hoseq = demo_hoseq()

    # phantom
    obj = brain_phantom2D(brain2D(); ss=simtype.ss, location=0.8); info(obj);
    obj.Δw .= simtype.B0 ? obj.Δw : obj.Δw * 0; # γ*1.5*(-3.45)*1e-6 * 2π
    obj.T2 .= simtype.T2 ? obj.T2 : obj.T2 * Inf; 

    BHO_recos = ["000", "100", "010", "001", "011", "101", "110", "111"]

    for idx in eachindex(BHO_recos)
        BHO_name = BHO_recos[idx]
        # scanner & sim_params
        sys = Scanner();
        sim_params = KomaMRICore.default_sim_params(); 
        sim_params["sim_method"] = BlochHighOrder(BHO_name);
        # sim_params["Nblocks"] = 150;
        sim_params["gpu"] = true;
        # sim_params["gpu_device"] = 1;
        sim_params["return_type"]="mat";

        # simulate
        signal = simulate(obj, hoseq, sys; sim_params);
        raw = signal_to_raw_data(signal, hoseq, :nominal)
        
        protocolName = "$(hoseq.SEQ.DEF["Name"])_$(BHO_name)"

        p = plot_image(reconstruct_2d_image(raw); title="$(sim_params["sim_method"])", height=700, width=750)
        savefig(p,  "$(path)/$(protocolName).svg",format="svg", height=700, width=750)
        raw.params["protocolName"] = protocolName
        mrd = ISMRMRDFile("$(path)/$(protocolName).mrd")
        save(mrd, raw)
    end
end


