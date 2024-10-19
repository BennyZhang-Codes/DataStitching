
"""
    raw = load_raw(BlochHighOrder("000"); simtype=SimType("B0wo_T2w_ss3"))
    load rawdata file from `src/demo/raw` folder and return a `RawAcquisitionData` object.

# Arguments
- `BHO`: (`::BlochHighOrder`) BlochHighOrder struct
# Keywords
- `simtype`: (`::SimType`) SimType struct

# Returns
- `raw`: (`::RawAcquisitionData`)
"""
function load_raw(BHO::BlochHighOrder; simtype::SimType=SimType("B0wo_T2w_ss3")) ::RawAcquisitionData
    folder = simtype.name
    raw_folder = "$(@__DIR__)/rawdata/$(folder)"
    @assert ispath(raw_folder) "folder not exist: $(raw_folder)"
    @info "rawdata for demo" sim_method=BHO.name mrd_file="$(@__DIR__)/raw/$(folder)/xw_sp2d-1mm-r1_$(BHO.name).mrd"

    raw_file = "$(@__DIR__)/rawdata/$(folder)/xw_sp2d-1mm-r1_$(BHO.name).mrd"
    @assert ispath(raw_file) "the raw file does not exist: $(raw_file)"
    raw = RawAcquisitionData(ISMRMRDFile(raw_file));
    return raw
end


"""
    seq = load_seq(; seq="spiral", r=2)
    load sequence file from `src/demo/seq` folder and return a `Sequence` object.
# Keywords
- `seqname`: (`::String`) Sequence name
- `r`: (`::Int64`) under sampling factor for Parallel Imaging

# Returns
- `seq`: (`::Sequence`)
"""
function load_seq(; seqname::String="demo", r::Int64=1)
    seq_path = "$(@__DIR__)/seq/$(seqname)_R$(r).seq"
    @assert ispath(seq_path) "the sequence file does not exist: $(seq_path)"
    seq = read_seq(seq_path)
    return seq
end

"""
    dfc_grad = load_dfc(; seq="spiral", r=2)
    load dfc file from `src/demo/dfc` folder and return a `Grad` object.
# Keywords
- `dfc_method`: (`::Symbol`) `:Stitched` or `:Standard`
- `seqname`: (`::String`) Sequence name
- `r`: (`::Int64`) under sampling factor for Parallel Imaging

# Returns
- `GR_dfc`: 
"""
function load_dfc(;dfc_method::Symbol=:Stitched, seqname::String="demo", r::Int64=1)
    dfc_path = "$(@__DIR__)/dfc/$(seqname)_R$(r).mat"
    @assert ispath(dfc_path) "the dfc file does not exist: $(dfc_path)"
    @assert dfc_method in [:Stitched, :Standard] "dfc_method must be :Stitched or :Standard"

    grad = MAT.matread(dfc_path);
    Δt = grad["dt"];
    nGradPoint, nTerm = size(grad["skopeStitched"])

    dfcStitched = [zeros(nTerm) grad["skopeStitched"]'] * 1e-3; 
    dfcStandard = [zeros(nTerm) grad["skopeStandard"]'] * 1e-3;
    ntStitched = grad["ntStitched"]
    ntStandard = grad["ntStandard"]

    t = Δt * ones(nGradPoint);
    if dfc_method == :Stitched
        GR_dfc = reshape([KomaMRIBase.Grad(dfcStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    elseif dfc_method == :Standard
        GR_dfc = reshape([KomaMRIBase.Grad(dfcStandard[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    end
    return GR_dfc, ntStitched, ntStandard
end

function load_dfc_mat(dfc_path::String)
    @assert ispath(dfc_path) "the dfc file does not exist: $(dfc_path)"
    grad = MAT.matread(dfc_path);
    Δt = grad["dt"];
    nGradPoint, nTerm = size(grad["skopeStitched"])

    dfcStitched = [zeros(nTerm) grad["skopeStitched"]'] * 1e-3; 
    dfcStandard = [zeros(nTerm) grad["skopeStandard"]'] * 1e-3;
    ntStitched = grad["ntStitched"]
    ntStandard = grad["ntStandard"]

    t = Δt * ones(nGradPoint);
    GR_dfcStitched = reshape([KomaMRIBase.Grad(dfcStitched[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    GR_dfcStandard = reshape([KomaMRIBase.Grad(dfcStandard[idx,:], t, 0, 0, 0) for idx=1:9], :, 1);
    return GR_dfcStitched, GR_dfcStandard, ntStitched, ntStandard
end
