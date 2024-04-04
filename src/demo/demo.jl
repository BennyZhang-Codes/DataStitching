export demo

function demo() ::HO_Sequence
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
    seq.GR[1,:] = -seq.GR[1,:]; 
    hoseq = HO_Sequence(seq);
    hoseq.GR_skope[2:4, :] = hoseq.SEQ.GR;
    hoseq.GR_skope[:,8] = GR_skope;
    return hoseq
end


