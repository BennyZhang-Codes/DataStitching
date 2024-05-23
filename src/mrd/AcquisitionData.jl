# import KomaHighOrder.KomaMRICore.KomaMRIBase: AcquisitionData
import KomaHighOrder.MRIBase: AcquisitionData, contrasts, slices, repetitions, trajectory, subsampleIndices, rawdata
"""
    AcquisitionData(f::RawAcquisitionData; estimateProfileCenter::Bool=false, OffsetBruker=false)

converts `RawAcquisitionData` into the equivalent `AcquisitionData` object.
If OffsetBruker=true, add a phase offset to correct position only along Y/Z direction
"""
function AcquisitionData(
    f::RawAcquisitionData, 
    sim_method::BlochHighOrder; 
    estimateProfileCenter::Bool=false, 
    sim_params=Dict{String,Any}())
    #Simulation parameter unpacking, and setting defaults if key is not defined
    sim_params = KomaMRICore.default_sim_params(sim_params)
    T = sim_params["precision"] == "f32" ? Float32 : Float64


    contrs = sort(unique(contrasts(f)))
    sls = sort(unique(slices(f)))
    reps = sort(unique(repetitions(f)))
    numContr = length(unique(contrasts(f)))
    numSl = length(unique(slices(f)))
    numRep = length(unique(repetitions(f)))

    # tr = [trajectory(f,contrast=contr) for contr=contrs]
    subsampleIdx = [subsampleIndices(f,contrast=contr,estimateProfileCenter=estimateProfileCenter) for contr=contrs]
    # kdata = [rawdata(f; contrast=i, slice=j, repetition=k) for i=contrs, j=sls, k=reps]

    tr = Trajectory{T}[]
    for contr=contrs
        t = trajectory(f,contrast=contr)
        t1 = Trajectory(T.(t.nodes), t.numProfiles, t.numSamplingPerProfile;
        times=T.(t.times), TE=T.(t.TE), AQ=T.(t.AQ), numSlices=t.numSlices, cartesian=t.cartesian, circular=t.circular)
        push!(tr, t1)
    end
    kdata = [Complex{T}.(rawdata(f; contrast=i, slice=j, repetition=k)) for i=contrs, j=sls, k=reps]
    acq = AcquisitionData(tr, kdata,
                            idx = subsampleIdx,
                            encodingSize = ntuple(d->f.params["encodedSize"][d], ndims(tr[1])),
                            fov = Float64.(ntuple(d->f.params["encodedFOV"][d], 3)))
    return acq
end

"""
t.AQ
t.TE
t.cartesian
t.circular
t.name
t.nodes
t.numProfiles
t.numSamplingPerProfile
t.numSlices
t.times

function Trajectory(nodes::AbstractMatrix{T}, numProfiles::Int64, numSamplingPerProfile;
    times=nothing, TE=0.0, AQ=1.e-3, numSlices::Int64=1,
    cartesian::Bool=false, circular::Bool=false) where T <: AbstractFloat

TE_ = T(TE)
AQ_ = T(AQ)

if times != nothing
ttimes = times
else
ttimes = readoutTimes(T,numProfiles,numSamplingPerProfile, numSlices; TE=TE_, AQ=AQ_)
end

return Trajectory("Custom", nodes, ttimes, TE_, AQ_, numProfiles, numSamplingPerProfile, numSlices, cartesian, circular)
end

"""