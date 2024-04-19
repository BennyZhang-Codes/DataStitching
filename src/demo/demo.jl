export demo, demo_obj, demo_seq, demo_GR_skope, demo_raw, demo_hoseq, demo_sim

function demo() ::Nothing
    @info "demos:"
    @info "    1. seq => demo_seq()"
    @info "    2. GR_skope => demo_GR_skope()"
    @info "    3. hoseq => demo_hoseq()"
    @info "    4. raw (mrd)  => demo_raw(\"000\")"
    @info "    5. sim (mat)  => raw, image = demo_sim()"
end

function demo_obj(; axis="axial", ss=5, location=0.8)
    @assert 0 <= location <= 1 "location must be between 0 and 1"
    path = "$(dirname(dirname(@__DIR__)))/$(phantom_dict[:path])/$(phantom_dict[:brain2d])"
    @assert isfile(path) "the phantom file does not exist: $(path)"
    data = MAT.matread(path)["data"]

    M, N, Z = size(data)
    if axis == "axial"
        z = Int32(ceil(Z*location))
        class = data[1:ss:end,1:ss:end, z]
    elseif axis == "coronal"
        m = Int32(ceil(M*location))
        class = data[m, 1:ss:end,1:ss:end]   
    elseif axis == "sagittal"
        n = Int32(ceil(N*location))
        class = data[1:ss:end, n,1:ss:end]
    end
    ρ = (class.==23)*1 .+ #CSF
        (class.==46)*.86 .+ #GM
        (class.==70)*.77 .+ #WM
        (class.==93)*1 .+ #FAT1
        (class.==116)*1 .+ #MUSCLE
        (class.==139)*.7 .+ #SKIN/MUSCLE
        (class.==162)*0 .+ #SKULL
        (class.==185)*0 .+ #VESSELS
        (class.==209)*.77 .+ #FAT2
        (class.==232)*1 .+ #DURA
        (class.==255)*.77 #MARROW
    return ρ
end

function demo_seq()
    path = @__DIR__
    seq = read_seq(path*"/xw_sp2d-1mm-r1_noDUM.seq");
    return seq
end

function demo_GR_skope()
    path = @__DIR__
    grad = MAT.matread(path*"/grad_1mm.mat");
    Δt = grad["dt"];
    skopeStitched = [zeros(9) grad["skopeStitched"]'] * 1e-3; 
    skopeStandard = [zeros(9) grad["skopeStandard"]'] * 1e-3;
    t = Δt * ones(88100);
    GR_skope = reshape([KomaMRIBase.Grad(skopeStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    return GR_skope
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

function demo_raw(name::String; folder::String="base") ::RawAcquisitionData
    @assert ispath("$(@__DIR__)/demo_raw/$(folder)") "folder not exist: $(folder)"
    BHO = BlochHighOrder(name)
    @info "demo_raw: $(BHO) $(@__DIR__)/demo_raw/$(folder)/xw_sp2d-1mm-r1_$(BHO.name)_nominal.mrd"

    raw_file = "$(@__DIR__)/demo_raw/$(folder)/xw_sp2d-1mm-r1_$(BHO.name)_nominal.mrd"
    @assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
    raw = RawAcquisitionData(ISMRMRDFile(raw_file));
    return raw
end


function demo_sim(;
    hoseq = demo_hoseq(),
    obj = brain_phantom2D(brain3D_02(); ss=3, location=0.8),
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