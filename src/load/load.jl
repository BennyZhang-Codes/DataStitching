export load_raw, load_seq, load_dfc
export load_dfc_mat
export load_hoseq


"""
    seq = load_seq(; seq="spiral", r=2)

# Description
    using read_seq function to load sequence file from `src/example/seq` folder and return a `Sequence` object.

# Keywords
- `seqname`: (`::String`) Sequence name
- `r`: (`::Int64`) under sampling factor for Parallel Imaging

# Returns
- `seq`: (`::Sequence`)
"""
function load_seq(; seqname::String="spiral", r::Int64=1)::Sequence
    seq_path = "$(@__DIR__)/seq/$(seqname)_R$(r).seq"
    @assert ispath(seq_path) "the sequence file does not exist: $(seq_path)"
    seq = read_seq(seq_path)
    return seq
end

"""
    GR_dfcStitched, GR_dfcStandard, ntStitched, ntStandard = load_dfc_mat(dfc_path)

# Description
    This function loads the dfc (dynamic field changes) file (`*.mat`) and returns the dynamic gradient of the stitched and standard methods.
the dfc file (*.mat) contains the following fields:
- `dt`               time interval between two time points [s].
- `ntStandard`       number of time points of the standard method.
- `ntStitched`       number of time points of the stitched method. A vector, one element for each segment.
- `skopeStandard`    dynamic measurment (up to 2nd order) of the standard method.
- `skopeStitched`    dynamic measurment (up to 2nd order) of the stitched method.

# Arguments
- `dfc_path`: (`::String`) path to the dfc file

# Returns
- `GR_dfcStitched`: (`::Array{Grad,2}`) stitched dynamic gradient
- `GR_dfcStandard`: (`::Array{Grad,2}`) standard dynamic gradient
- `ntStitched`: (`::Int64`) number of time points of the stitched method
- `ntStandard`: (`::Int64`) number of time points of the standard method
"""
function load_dfc_mat(dfc_path::String)
    @assert ispath(dfc_path) "the dfc file does not exist: $(dfc_path)"
    grad = MAT.matread(dfc_path);
    Δt = grad["dt"];
    nGradPoint, nTerm = size(grad["skopeStitched"])              
    dfcStitched = grad["skopeStitched"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
    dfcStandard = grad["skopeStandard"]' * 1e-3; # mT, mT/m, mT/m² => T, T/m, T/m²
    ntStitched  = grad["ntStitched"]
    ntStandard  = grad["ntStandard"]
    t = Δt * (nGradPoint-1);
    GR_dfcStitched = reshape([KomaMRIBase.Grad(dfcStitched[idx,:], t, Δt/2, Δt/2, 0) for idx=1:9], :, 1);
    GR_dfcStandard = reshape([KomaMRIBase.Grad(dfcStandard[idx,:], t, Δt/2, Δt/2, 0) for idx=1:9], :, 1);
    return GR_dfcStitched, GR_dfcStandard, ntStitched, ntStandard
end


"""
    dfc_grad = load_dfc(; seq="spiral", r=2)

# Description
    load dfc file from `src/example/dfc` folder and return a `Grad` object.

# Keywords
- `dfc_method`: (`::Symbol`) `:Stitched` or `:Standard`
- `seqname`: (`::String`) prefix of the sequence name: "spiral" or "gre" ...
- `r`: (`::Int64`) under sampling factor for Parallel Imaging

# Returns
- `GR_dfc`: 
"""
function load_dfc(;dfc_method::Symbol=:Stitched, seqname::String="spiral", r::Int64=1)
    dfc_path = "$(@__DIR__)/dfc/$(seqname)_R$(r).mat"
    @assert dfc_method in [:Stitched, :Standard] "dfc_method must be :Stitched or :Standard"

    GR_dfcStitched, GR_dfcStandard, ntStitched, ntStandard = load_dfc_mat(dfc_path);

    if dfc_method == :Stitched
        GR_dfc = GR_dfcStitched;
    elseif dfc_method == :Standard
        GR_dfc = GR_dfcStandard;
    end
    return GR_dfc, ntStitched, ntStandard
end

"""
    hoseq = load_hoseq(;dfc_method=:Stitched, seqname="spiral", r=2)

# Description
    This function loads the sequence (*.seq) and dfc file (*.mat) and returns a `HO_Sequence` object.

# Keywords
- `dfc`: (`::Bool`) whether to include dfc in the sequence
- `dfc_method`: (`::Symbol`) `:Stitched` or `:Standard`
- `seqname`: (`::String`) Sequence name
- `r`: (`::Int64`) under sampling factor for Parallel Imaging

# Returns
- `hoseq`: (`::HO_Sequence`) the HO_Sequence object
"""
function load_hoseq(;dfc_method::Symbol=:Stitched, seqname::String="spiral", r::Int64=1) ::HO_Sequence
    seq = load_seq(;seqname=seqname, r=r); # dfc sequence
    seq.GR[1,:] = -seq.GR[1,:];            # reverse the sign of the gradient (axis x)
    hoseq = HO_Sequence(seq);              # convert to HO_Sequence object
    GR_dfc, ntStitched, ntStandard = load_dfc(;dfc_method=dfc_method, seqname=seqname, r=r);
    hoseq.GR_dfc[2:4, :] = hoseq.SEQ.GR;
    hoseq.GR_dfc[:,8] = GR_dfc;            # "8" is the index of the readout block in the spiral sequence
    return hoseq
end
