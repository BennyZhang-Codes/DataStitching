export demo, demo_raw, demo_hoseq, demo_sim

function demo() ::Nothing
    @info "demos:"
    @info "    1. demo_hoseq => demo_hoseq()"
    @info "    2. raw (mrd)  => demo_raw(\"000\")"
end

function demo_hoseq() ::HO_Sequence
    path = @__DIR__
    seq = read_seq(path*"/xw_sp2d-1mm-r1_noDUM.seq"); # skope sequence
    grad = MAT.matread(path*"/grad_1mm.mat"); # skope measured gradients
    
    # skope 
    Δt = grad["dt"];
    skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3; 
    skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3;
    t = Δt * ones(88100);
    GR_skope = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    
    # hoseq
    seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
    hoseq = HO_Sequence(seq);
    hoseq.GR_skope[2:4, :] = hoseq.SEQ.GR;
    hoseq.GR_skope[:,8] = GR_skope;
    return hoseq
end

function demo_raw(name::String) ::RawAcquisitionData
    @assert name == "000" || name == "100" || name == "010" || name == "001" || name == "110" || name == "101" || 
    name == "011" || name == "111" "has raw: 000, 100, 010, 001, 110, 101, 011, 111"
    path = @__DIR__
    raw_file = path*"/demo_raw/xw_sp2d-1mm-r1_$(name)_nominal.mrd"
    @assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
    raw = RawAcquisitionData(ISMRMRDFile(raw_file));
    return raw
end


function demo_sim(;
    hoseq = demo_hoseq(),
    obj = brain_phantom2D(brain3D_02(); ss=3, location=0.8),
    sim_params=Dict{String,Any}(),
    sim_method::BlochHighOrder=BlochHighOrder("000")) ::Nothing
    path = @__DIR__
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
    return raw, reconstruct_2d_image(raw)
end