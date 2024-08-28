function get_kdata(raw::RawAcquisitionData, shape::Tuple)
    profiles = raw.profiles
    data = Array{Complex{Float32},length(shape)}(undef,shape);
    for idx in eachindex(profiles)
        profile = profiles[idx]
        header = profile.head
        # Stuff into the buffer
        Cha = Int(profile.head.active_channels)
        Z   = 1 + header.idx.kspace_encode_step_2
        Y   = 1 + header.idx.kspace_encode_step_1
        X   = Int(profile.head.number_of_samples)
        Avg = 1 + header.idx.average
        Sli = 1 + header.idx.slice
        Con = 1 + header.idx.contrast
        Pha = 1 + header.idx.phase
        Rep = 1 + header.idx.repetition
        Set = 1 + header.idx.set
        Seg = 1 + header.idx.segment
        data[:, Z, Y, :, Avg, Sli, Con, Pha, Rep, Set, Seg] = transpose(profile.data);
    end
    return data
end

# kdata = get_kdata(raw, shape)
# kdata = dropdims(kdata, dims = tuple(findall(size(kdata) .== 1)...))
# kdims = [dims[idx] for idx in 1:length(shape) if shape[idx]>1]