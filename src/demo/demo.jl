export demo, demo_seq, demo_GR_skope, demo_raw, demo_hoseq, demo_sim

function demo() ::Nothing
    @info "demos:"
    @info "    1. seq => demo_seq()"
    @info "    2. GR_skope => demo_GR_skope()"
    @info "    3. hoseq => demo_hoseq()"
    @info "    4. raw (mrd)  => demo_raw(\"000\")"
    @info "    5. sim (mat)  => raw, image = demo_sim()"
end

Base.@kwdef struct SimType
    B0::Bool = true
    T2::Bool = true
    ss::Int64 = 3
    name::String = "B0$(B0 ? "w" : "wo")_T2$(T2 ? "w" : "wo")_ss$(ss)"
end
export SimType
Base.show(io::IO, b::SimType) = begin
	print(io, "SimType[[$(b.name)] B0=$(b.B0) | T2=$(b.T2) | ss=$(b.ss) ]")
end

function SimType(name::String)
    @assert length(split(name,"_")) == 3 "Invalid name for SimType. Valid name like: B0wo_T2w_ss3"
    B0, T2, ss = split(name,"_")
    @assert !isnothing(findfirst("B0", B0)) && !isnothing(findfirst("T2", T2)) && !isnothing(findfirst("ss", ss)) "Invalid name for SimType. Valid name like: B0wo_T2w_ss3"
    B0, T2, ss = split(B0, "B0")[2], split(T2, "T2")[2], split(ss, "ss")[2]
    B0, T2, ss = B0 == "w", T2 == "w", parse(Int, ss)
    return SimType(B0=B0, T2=T2, ss=ss)
end

include("demo_raw/generate_raw.jl")
export generate_raw

function demo_seq()
    path = @__DIR__
    seq = read_seq(path*"/xw_sp2d-1mm-r1_noDUM.seq");
    return seq
end

function demo_GR_skope(;key::Symbol=:Stitched)
    @assert key in [:Stitched, :Standard] "key must be :Stitched or :Standard"
    path = @__DIR__
    grad = MAT.matread(path*"/grad_1mm.mat");
    Δt = grad["dt"];
    skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3; 
    skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3;
    t = Δt * ones(88100);
    if key == :Stitched
        GR_skope = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    elseif key == :Standard
        GR_skope = reshape([KomaMRIBase.Grad(skopeStandard[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    end
    return GR_skope
end


function demo_hoseq(;skope::Bool=true, skope_method::Symbol=:Stitched) ::HO_Sequence
    seq = demo_seq(); # skope sequence
    seq.GR[1,:] = -seq.GR[1,:]; # reverse the sign of the gradient (axis x)
    hoseq = HO_Sequence(seq); # hoseq
    if skope
        GR_skope = demo_GR_skope(;key=skope_method);
        hoseq.GR_skope[2:4, :] = hoseq.SEQ.GR;
        hoseq.GR_skope[:,8] = GR_skope;
    end
    return hoseq
end

function demo_raw(BHO::BlochHighOrder; simtype::SimType=SimType("B0wo_T2w_ss3")) ::RawAcquisitionData
    folder = simtype.name
    @assert ispath("$(@__DIR__)/demo_raw/$(folder)") "folder not exist: $(folder)"
    @info "demo_raw" sim_method=BHO.name mrd_file="$(@__DIR__)/demo_raw/$(folder)/xw_sp2d-1mm-r1_$(BHO.name).mrd"

    raw_file = "$(@__DIR__)/demo_raw/$(folder)/xw_sp2d-1mm-r1_$(BHO.name).mrd"
    @assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
    raw = RawAcquisitionData(ISMRMRDFile(raw_file));
    return raw
end


function demo_sim(;
    hoseq = demo_hoseq(),
    obj = brain_phantom2D(brain2D(); ss=3, location=0.8),
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