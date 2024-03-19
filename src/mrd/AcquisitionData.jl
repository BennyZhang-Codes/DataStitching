# import KomaHighOrder.KomaMRICore.KomaMRIBase: AcquisitionData
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory
"""
    AcquisitionData(f::RawAcquisitionData; estimateProfileCenter::Bool=false, OffsetBruker=false)

converts `RawAcquisitionData` into the equivalent `AcquisitionData` object.
If OffsetBruker=true, add a phase offset to correct position only along Y/Z direction
"""
function AcquisitionData(f::RawAcquisitionData, sim_method::BlochHighOrder; estimateProfileCenter::Bool=false)
    contrs = sort(unique(contrasts(f)))
    sls = sort(unique(slices(f)))
    reps = sort(unique(repetitions(f)))
    numContr = length(unique(contrasts(f)))
    numSl = length(unique(slices(f)))
    numRep = length(unique(repetitions(f)))
    # tr = [trajectory(f,contrast=contr) for contr=1:numContr]
    # subsampleIdx = [subsampleIndices(f,contrast=contr,estimateProfileCenter=estimateProfileCenter) for contr=1:numContr]
    # kdata = [rawdata(f; contrast=i, slice=j, repetition=k) for i=1:numContr, j=1:numSl, k=1:numRep]
    tr = [trajectory(f,contrast=contr) for contr=contrs]
    subsampleIdx = [subsampleIndices(f,contrast=contr,estimateProfileCenter=estimateProfileCenter) for contr=contrs]
    kdata = [rawdata(f; contrast=i, slice=j, repetition=k) for i=contrs, j=sls, k=reps]
    acq = AcquisitionData(tr, kdata,
                            idx = subsampleIdx,
                            encodingSize = ntuple(d->f.params["encodedSize"][d], ndims(tr[1])),
                            fov = Float64.(ntuple(d->f.params["encodedFOV"][d], 3)))
    return acq
end