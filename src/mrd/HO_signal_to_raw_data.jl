import KomaMRI.KomaMRICore: signal_to_raw_data
"""
    raw = signal_to_raw_data(signal, hoseq, sim_method; phantom_name, sys, sim_params)

Transforms the raw signal into a RawAcquisitionData struct (nearly equivalent to the ISMRMRD
format) used for reconstruction with MRIReco.

# Arguments
- `signal`: (`::Matrix{Complex}`) raw signal matrix
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct
- `sim_method`: (`::BlochHighOrder`) simulation method

# Keywords
- `phantom_name`: (`::String`, `="Phantom"`) phantom name
- `sys`: (`::Scanner`, `=Scanner()`) Scanner struct
- `sim_params`: (`::Dict{String, Any}`, `=Dict{String,Any}()`) simulation parameter dictionary

# Returns
- `raw`: (`::RawAcquisitionData`) RawAcquisitionData struct
"""
function signal_to_raw_data(
    signal, hoseq::HO_Sequence, sim_method::BlochHighOrder;
    phantom_name="Phantom", sys=Scanner(), sim_params=Dict{String,Any}(), ndims=6
)
    seq = hoseq.SEQ
    version = string(VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["version"]))
    #Number of samples and FOV
    _, ktraj, _, ktraj_dfc = get_kspace(hoseq) #kspace information
    mink = minimum(ktraj, dims=1)
    maxk = maximum(ktraj, dims=1)
    Wk = maxk .- mink
    Δx = 1 ./ Wk[1:2] #[m] Only x-y
    Nx = get(seq.DEF, "Nx", 1)
    Ny = get(seq.DEF, "Ny", 1)
    Nz = get(seq.DEF, "Nz", 1)

    adc_length        = UInt32(get(seq.DEF, "adc_length", maximum(seq.ADC.N)))
    adc_segmentLength = UInt32(get(seq.DEF, "adc_segmentLength", 1000))
    adc_nSegments     = UInt32(get(seq.DEF, "adc_nSegments", Int(ceil(adc_length/adc_segmentLength))))
    


    if haskey(seq.DEF, "FOV")
        FOVx, FOVy, _ = seq.DEF["FOV"] #[m]
        if FOVx > 1 FOVx *= 1e-3 end #mm to m, older versions of Pulseq saved FOV in mm
        if FOVy > 1 FOVy *= 1e-3 end #mm to m, older versions of Pulseq saved FOV in mm
        Nx = round(Int64, FOVx / Δx[1])
        Ny = round(Int64, FOVy / Δx[2])
    else
        FOVx = Nx * Δx[1]
        FOVy = Ny * Δx[2]
    end
    #It needs to be transposed for the raw data
    ktraj = maximum(2*abs.(ktraj[:])) == 0 ? transpose(ktraj) : transpose(ktraj)./ maximum(2*abs.(ktraj[:]))
    ktraj_dfc = maximum(2*abs.(ktraj_dfc[:])) == 0 ? transpose(ktraj_dfc) : transpose(ktraj_dfc)./ maximum(2*abs.(ktraj_dfc[:]))
    ktraj = [ktraj; ktraj_dfc]

    #First we define the ISMRMRD data XML header
    #userParameters <- sim_params
    for (key, val) in sim_params
        if typeof(val) <: Integer #Fixes problem with bools
            sim_params[key] = Int(val)
        end
    end
    #XML header
    params = Dict(
        #AcquisitionSystemInformation
        "systemVendor"                   => "KomaHighOrder.jl", #String
        "systemModel"                    => "v"*version, #String
        "systemFieldStrength_T"          => sys.B0, #Float
        "institutionName"                => "BennyZhang", #String
        #subjectInformation
        "patientName"                    => phantom_name,
        #experimentalConditions
        "H1resonanceFrequency_Hz"        => floor(Int64, γ * sys.B0), #Long (Int)
        #measurementInformation
        "protocolName"                   => haskey(seq.DEF,"Name") ? seq.DEF["Name"] : "NoName", #String
        # "trajectoryDescription"          => Dict{String, Any}("comment"=>""), #You can put wathever you want here: comment, bandwidth, MaxGradient_G_per_cm, MaxSlewRate_G_per_cm_per_s, interleaves, etc
        #encoding
        #   encodedSpace
        "encodedSize"                    => [Nx, Ny, 1],                        #encodedSpace>matrixSize
        "encodedFOV"                     => Float32.([FOVx, FOVy, 1e-3]*1e3),   #encodedSpace>fieldOfView_mm
        #   reconSpace
        "reconSize"                      => [Nx+Nx%2, Ny+Ny%2, 1],              #reconSpace>matrixSize
        "reconFOV"                       => Float32.([FOVx, FOVy, 1e-3]*1e3),   #reconSpace>fieldOfView_mm
        #encodingLimits
        "enc_lim_kspace_encoding_step_0" => Limit(0, Nx-1, ceil(Int, Nx / 2)),  #min, max, center, e.g. phase encoding line number
        "enc_lim_kspace_encoding_step_1" => Limit(0, Ny-1, ceil(Int, Ny / 2)),  #min, max, center, e.g. partition encoding number
        "enc_lim_kspace_encoding_step_2" => Limit(0, 0, 0),                     #min, max, center, e.g. partition encoding number
        "enc_lim_average"                => Limit(0, 0, 0),                     #min, max, center, e.g. signal average number
        "enc_lim_slice"                  => Limit(0, 0, 0),                     #min, max, center, e.g. imaging slice number
        "enc_lim_contrast"               => Limit(0, 0, 0),                     #min, max, center, e.g. echo number in multi-echo
        "enc_lim_phase"                  => Limit(0, 0, 0),                     #min, max, center, e.g. cardiac phase number
        "enc_lim_repetition"             => Limit(0, 0, 0),                     #min, max, center, e.g. dynamic number for dynamic scanning
        "enc_lim_set"                    => Limit(0, 0, 0),                     #min, max, center, e.g. flow encoding set
        "enc_lim_segment"                => Limit(0, 0, 0),                     #min, max, center, e.g. segment number for segmented acquisition
        "trajectory"                     => "other",
        #sequenceParameters
        # "TR"                             => 0,
        # "TE"                             => 0,
        # "TI"                             => 0,
        # "flipAngle_deg"                  => 0,
        # "echo_spacing"                   => 0,
        "userParameters"                 => sim_params, #Dict with parameters
    )

    #Then, we define the Profiles
    profiles = Profile[]
    t_acq = get_adc_sampling_times(seq)
    Nadcs = sum(is_ADC_on.(seq))
    NadcsPerImage = floor(Int, Nadcs / Nz)
    scan_counter = 0
    nz = 0
    current = 1


    # adc_length        = get(seq.DEF, "adc_length", maximum(seq.ADC.N))
    # adc_segmentLength = get(seq.DEF, "adc_segmentLength", 1000)
    # adc_nSegments     = get(seq.DEF, "adc_nSegments", Int(ceil(adc_length/adc_segmentLength)))

    for s = seq #Iterate over sequence blocks
        if is_ADC_on(s)
            Nsamples = s.ADC.N[1]
            Δt_us = floor( s.ADC.T[1] / (Nsamples - 1) * 1e6 )
            # t0_us = floor( t_acq[current] * 1e6 )
            flag  = 0
            if scan_counter == 0
                flag += KomaMRICore.ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1
                flag += KomaMRICore.ISMRMRD_ACQ_FIRST_IN_SLICE
            elseif scan_counter == Nadcs - 1
                flag += KomaMRICore.ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1
                flag += KomaMRICore.ISMRMRD_ACQ_LAST_IN_SLICE
            end

            for seg = 0:adc_nSegments-1
                sample_range_start = current + seg * adc_segmentLength
                sample_range_end = adc_nSegments != 1 ? min(adc_length, current + (seg + 1) * adc_segmentLength - 1) : current + seg * adc_segmentLength + Nsamples - 1
                number_of_samples = sample_range_end - sample_range_start + 1

                t0_us = floor( t_acq[sample_range_start] * 1e6 )
                #Header of profile data, head::AcquisitionHeader
                head = AcquisitionHeader(
                    UInt16(1), #version uint16: First unsigned int indicates the version
                    UInt64(flag), #flags uint64: bit field with flags, noise measure, calibration, coil sens, measure, etc
                    UInt32(0), #measurement_uid uint32: Unique ID for the measurement
                    UInt32(scan_counter+seg), #scan_counter uint32: Current acquisition number in the measurement
                    UInt32(t0_us), #acquisition_time_stamp uint32: Acquisition clock, I am "miss"-using this variable to store t0 in us
                    UInt32.((0, 0, 0)), #physiology_time_stamp uint32x3: Physiology time stamps, e.g. ecg, breating, etc.
                    UInt16(number_of_samples), #number_of_samples uint16
                    UInt16(1), #available_channels uint16: Available coils
                    UInt16(1), #active_channels uint16: Active coils on current acquisiton
                    Tuple(UInt64(0) for i=1:16), #channel_mask uint64x16: Active coils on current acquisiton
                    UInt16(0), #discard_pre uint16: Samples to be discarded at the beginning of acquisition
                    UInt16(0), #discard_post uint16: Samples to be discarded at the end of acquisition
                    UInt16(0), #center_sample uint16: Sample at the center of k-space
                    UInt16(0), #encoding_space_ref uint16: Reference to an encoding space, typically only one per acquisition
                    UInt16(ndims), #trajectory_dimensions uint16: Indicates the dimensionality of the trajectory vector (0 means no trajectory)
                    Float32(Δt_us), #sample_time_us float32: Time between samples in micro seconds, sampling BW
                    Float32.((0, 0, 0)), #position float32x3: Three-dimensional spatial offsets from isocenter
                    Float32.((1, 0, 0)), #read_dir float32x3: Directional cosines of the readout/frequency encoding
                    Float32.((0, 1, 0)), #phase_dir float32x3: Directional cosines of the phase
                    Float32.((0, 0, 1)), #slice_dir float32x3: Directional cosines of the slice direction
                    Float32.((0, 0, 0)), #patient_table_position float32x3: Patient table off-center
                    EncodingCounters( #idx uint16x17: Encoding loop counters
                        UInt16(scan_counter+seg), #kspace_encode_step_1 uint16: e.g. phase encoding line number
                        UInt16(0), #kspace_encode_step_2 uint16: e.g. partition encoding number
                        UInt16(0), #average uint16: e.g. signal average number
                        UInt16(nz), #slice uint16: e.g. imaging slice number
                        UInt16(0), #contrast uint16: e.g. echo number in multi-echo
                        UInt16(0), #phase uint16: e.g. cardiac phase number
                        UInt16(0), #repetition uint16: e.g. dynamic number for dynamic scanning
                        UInt16(0), #set uint16: e.g. flow encoding set
                        UInt16(seg), #segment uint16: e.g. segment number for segmented acquisition
                        Tuple(UInt16(0) for i=1:8) #user uint16x8: Free user parameters
                    ),
                    Tuple(Int32(0) for i=1:8), #user_int int32x8: Free user parameters
                    Tuple(Float32(0) for i=1:8) #user_float float32x8: Free user parameters
                )
                # #Trajectory information, traj::Array{Float32,2}, 1dim=DIM, 2dim=numsaples
                # traj = ktraj[1:ndims, current:current+Nsamples-1]
                # #Acquired data, data::Array{Complex{Float32},2}, 1dim=numsamples, 2dim=coils
                # dat =  signal[current:current+Nsamples-1, :]
                #Trajectory information, traj::Array{Float32,2}, 1dim=DIM, 2dim=numsaples

                traj = ktraj[1:ndims, sample_range_start:sample_range_end]
                #Acquired data, data::Array{Complex{Float32},2}, 1dim=numsamples, 2dim=coils
                dat =  signal[sample_range_start:sample_range_end, :]

                #Saving profile
                push!(profiles, Profile(head, Float32.(traj), ComplexF32.(dat)))
            end
            #Update counters
            # scan_counter += 1
            scan_counter += adc_nSegments

            current += Nsamples
            if scan_counter % NadcsPerImage == 0 #For now only Nz is considered
                nz += 1 #another image
                scan_counter = 0 #reset counter
            end
        end
    end

    return RawAcquisitionData(params, profiles)
end


"""
    raw = signal_to_raw_data(signal, hoseq, key; phantom_name, sys, sim_params)

Transforms the raw signal into a RawAcquisitionData struct (nearly equivalent to the ISMRMRD
format) used for reconstruction with MRIReco.

# Arguments
- `signal`: (`::Matrix{Complex}`) raw signal matrix
- `hoseq`: (`::HO_Sequence`) HO_Sequence struct
- `key`: (`::Symbol`) simulation method

# Keywords
- `phantom_name`: (`::String`, `="Phantom"`) phantom name
- `sys`: (`::Scanner`, `=Scanner()`) Scanner struct
- `sim_params`: (`::Dict{String, Any}`, `=Dict{String,Any}()`) simulation parameter dictionary

# Returns
- `raw`: (`::RawAcquisitionData`) RawAcquisitionData struct
"""
function signal_to_raw_data(
    signal::Any, hoseq::HO_Sequence, key::Symbol;
    phantom_name="Phantom", sys=Scanner(), sim_params=Dict{String,Any}(), ndims=2
)
    if size(signal, 3) == 1
        signal = signal[:,:,1]
    end
    return signal_to_raw_data(signal, hoseq, key, phantom_name=phantom_name, sys=sys, sim_params=sim_params, ndims=ndims)
end
function signal_to_raw_data(
    signal::Matrix, hoseq::HO_Sequence, key::Symbol;
    phantom_name="Phantom", sys=Scanner(), sim_params=Dict{String,Any}(), ndims=2
)
    sim_params = KomaMRICore.default_sim_params(sim_params)
    TotalSamples, Ncoils = size(signal)
    @assert key == :nominal || key == :measured  "key must be :nominal or :measured"
    seq = hoseq.SEQ
    version = string(VersionNumber(Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["version"]))
    #Number of samples and FOV
    _, ktraj, _, ktraj_dfc = get_kspace(hoseq) #kspace information #kspace information
    if key == :nominal
        ktraj = ktraj
    elseif key == :measured
        ktraj = ktraj_dfc[:,2:4]
    end

    if sim_params["precision"] == "f32" #Default
        signal = signal |> f32 
        ktraj  = ktraj  |> f32 
    elseif sim_params["precision"] == "f64"
        signal = signal |> f64
        ktraj  = ktraj  |> f64
    end
    mink = minimum(ktraj, dims=1)
    maxk = maximum(ktraj, dims=1)
    Wk = maxk .- mink
    Δx = 1 ./ Wk[1:2] #[m] Only x-y
    Nx = get(seq.DEF, "Nx", 1)
    Ny = get(seq.DEF, "Ny", 1)
    Nz = get(seq.DEF, "Nz", 1)

    adc_length        = UInt32(get(seq.DEF, "adc_length", maximum(seq.ADC.N)))
    adc_segmentLength = UInt32(get(seq.DEF, "adc_segmentLength", 1000))
    adc_nSegments     = UInt32(get(seq.DEF, "adc_nSegments", Int(ceil(adc_length/adc_segmentLength))))
    


    if haskey(seq.DEF, "FOV")
        FOVx, FOVy, _ = seq.DEF["FOV"] #[m]
        if FOVx > 1 FOVx *= 1e-3 end #mm to m, older versions of Pulseq saved FOV in mm
        if FOVy > 1 FOVy *= 1e-3 end #mm to m, older versions of Pulseq saved FOV in mm
        Nx = round(Int64, FOVx / Δx[1])
        Ny = round(Int64, FOVy / Δx[2])
    else
        FOVx = Nx * Δx[1]
        FOVy = Ny * Δx[2]
    end
    #It needs to be transposed for the raw data
    ktraj = maximum(2*abs.(ktraj[:])) == 0 ? transpose(ktraj) : transpose(ktraj)./ maximum(2*abs.(ktraj[:]))

    #First we define the ISMRMRD data XML header
    #userParameters <- sim_params
    for (key, val) in sim_params
        if typeof(val) <: Integer #Fixes problem with bools
            sim_params[key] = Int(val)
        end
    end
    #XML header
    params = Dict(
        #AcquisitionSystemInformation
        "systemVendor"                   => "KomaHighOrder.jl", #String
        "systemModel"                    => "v"*version, #String
        "systemFieldStrength_T"          => sys.B0, #Float
        "institutionName"                => "BennyZhang", #String
        #subjectInformation
        "patientName"                    => phantom_name,
        #experimentalConditions
        "H1resonanceFrequency_Hz"        => floor(Int64, γ * sys.B0), #Long (Int)
        #measurementInformation
        "protocolName"                   => haskey(seq.DEF,"Name") ? seq.DEF["Name"] : "NoName", #String
        # "trajectoryDescription"          => Dict{String, Any}("comment"=>""), #You can put wathever you want here: comment, bandwidth, MaxGradient_G_per_cm, MaxSlewRate_G_per_cm_per_s, interleaves, etc
        #encoding
        #   encodedSpace
        "encodedSize"                    => [Nx, Ny, 1],                        #encodedSpace>matrixSize
        "encodedFOV"                     => Float32.([FOVx, FOVy, 1e-3]*1e3),   #encodedSpace>fieldOfView_mm
        #   reconSpace
        "reconSize"                      => [Nx+Nx%2, Ny+Ny%2, 1],              #reconSpace>matrixSize
        "reconFOV"                       => Float32.([FOVx, FOVy, 1e-3]*1e3),   #reconSpace>fieldOfView_mm
        #encodingLimits
        "enc_lim_kspace_encoding_step_0" => Limit(0, Nx-1, ceil(Int, Nx / 2)),  #min, max, center, e.g. phase encoding line number
        "enc_lim_kspace_encoding_step_1" => Limit(0, Ny-1, ceil(Int, Ny / 2)),  #min, max, center, e.g. partition encoding number
        "enc_lim_kspace_encoding_step_2" => Limit(0, 0, 0),                     #min, max, center, e.g. partition encoding number
        "enc_lim_average"                => Limit(0, 0, 0),                     #min, max, center, e.g. signal average number
        "enc_lim_slice"                  => Limit(0, 0, 0),                     #min, max, center, e.g. imaging slice number
        "enc_lim_contrast"               => Limit(0, 0, 0),                     #min, max, center, e.g. echo number in multi-echo
        "enc_lim_phase"                  => Limit(0, 0, 0),                     #min, max, center, e.g. cardiac phase number
        "enc_lim_repetition"             => Limit(0, 0, 0),                     #min, max, center, e.g. dynamic number for dynamic scanning
        "enc_lim_set"                    => Limit(0, 0, 0),                     #min, max, center, e.g. flow encoding set
        "enc_lim_segment"                => Limit(0, 0, 0),                     #min, max, center, e.g. segment number for segmented acquisition
        "trajectory"                     => "other",
        #sequenceParameters
        # "TR"                             => 0,
        # "TE"                             => 0,
        # "TI"                             => 0,
        # "flipAngle_deg"                  => 0,
        # "echo_spacing"                   => 0,
        "userParameters"                 => sim_params, #Dict with parameters
    )

    #Then, we define the Profiles
    profiles = Profile[]
    t_acq = get_adc_sampling_times(seq)
    Nadcs = sum(is_ADC_on.(seq))
    NadcsPerImage = floor(Int, Nadcs / Nz)
    scan_counter = 0
    nz = 0
    current = 1


    # adc_length        = get(seq.DEF, "adc_length", maximum(seq.ADC.N))
    # adc_segmentLength = get(seq.DEF, "adc_segmentLength", 1000)
    # adc_nSegments     = get(seq.DEF, "adc_nSegments", Int(ceil(adc_length/adc_segmentLength)))

    for s = seq #Iterate over sequence blocks
        if is_ADC_on(s)
            Nsamples = s.ADC.N[1]
            Δt_us = floor( s.ADC.T[1] / (Nsamples - 1) * 1e6 )
            # t0_us = floor( t_acq[current] * 1e6 )
            flag  = 0
            if scan_counter == 0
                flag += KomaMRICore.ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1
                flag += KomaMRICore.ISMRMRD_ACQ_FIRST_IN_SLICE
            elseif scan_counter == Nadcs - 1
                flag += KomaMRICore.ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1
                flag += KomaMRICore.ISMRMRD_ACQ_LAST_IN_SLICE
            end

            for seg = 0:adc_nSegments-1
                sample_range_start = current + seg * adc_segmentLength
                # sample_range_end = adc_nSegments != 1 ? min(adc_length, current + (seg + 1) * adc_segmentLength - 1) : current + seg * adc_segmentLength + Nsamples - 1
                sample_range_end = sample_range_start +  adc_segmentLength - 1
                number_of_samples = sample_range_end - sample_range_start + 1

                t0_us = floor( t_acq[sample_range_start] * 1e6 )
                #Header of profile data, head::AcquisitionHeader
                head = AcquisitionHeader(
                    UInt16(1), #version uint16: First unsigned int indicates the version
                    UInt64(flag), #flags uint64: bit field with flags, noise measure, calibration, coil sens, measure, etc
                    UInt32(0), #measurement_uid uint32: Unique ID for the measurement
                    UInt32(scan_counter+seg), #scan_counter uint32: Current acquisition number in the measurement
                    UInt32(t0_us), #acquisition_time_stamp uint32: Acquisition clock, I am "miss"-using this variable to store t0 in us
                    UInt32.((0, 0, 0)), #physiology_time_stamp uint32x3: Physiology time stamps, e.g. ecg, breating, etc.
                    UInt16(number_of_samples), #number_of_samples uint16
                    UInt16(Ncoils), #available_channels uint16: Available coils
                    UInt16(Ncoils), #active_channels uint16: Active coils on current acquisiton
                    Tuple(UInt64(0) for i=1:16), #channel_mask uint64x16: Active coils on current acquisiton
                    UInt16(0), #discard_pre uint16: Samples to be discarded at the beginning of acquisition
                    UInt16(0), #discard_post uint16: Samples to be discarded at the end of acquisition
                    UInt16(0), #center_sample uint16: Sample at the center of k-space
                    UInt16(0), #encoding_space_ref uint16: Reference to an encoding space, typically only one per acquisition
                    UInt16(ndims), #trajectory_dimensions uint16: Indicates the dimensionality of the trajectory vector (0 means no trajectory)
                    Float32(Δt_us), #sample_time_us float32: Time between samples in micro seconds, sampling BW
                    Float32.((0, 0, 0)), #position float32x3: Three-dimensional spatial offsets from isocenter
                    Float32.((1, 0, 0)), #read_dir float32x3: Directional cosines of the readout/frequency encoding
                    Float32.((0, 1, 0)), #phase_dir float32x3: Directional cosines of the phase
                    Float32.((0, 0, 1)), #slice_dir float32x3: Directional cosines of the slice direction
                    Float32.((0, 0, 0)), #patient_table_position float32x3: Patient table off-center
                    EncodingCounters( #idx uint16x17: Encoding loop counters
                        UInt16(scan_counter+seg), #kspace_encode_step_1 uint16: e.g. phase encoding line number
                        UInt16(0), #kspace_encode_step_2 uint16: e.g. partition encoding number
                        UInt16(0), #average uint16: e.g. signal average number
                        UInt16(nz), #slice uint16: e.g. imaging slice number
                        UInt16(0), #contrast uint16: e.g. echo number in multi-echo
                        UInt16(0), #phase uint16: e.g. cardiac phase number
                        UInt16(0), #repetition uint16: e.g. dynamic number for dynamic scanning
                        UInt16(0), #set uint16: e.g. flow encoding set
                        UInt16(seg), #segment uint16: e.g. segment number for segmented acquisition
                        Tuple(UInt16(0) for i=1:8) #user uint16x8: Free user parameters
                    ),
                    Tuple(Int32(0) for i=1:8), #user_int int32x8: Free user parameters
                    Tuple(Float32(0) for i=1:8) #user_float float32x8: Free user parameters
                )
                # #Trajectory information, traj::Array{Float32,2}, 1dim=DIM, 2dim=numsaples
                # traj = ktraj[1:ndims, current:current+Nsamples-1]
                # #Acquired data, data::Array{Complex{Float32},2}, 1dim=numsamples, 2dim=coils
                # dat =  signal[current:current+Nsamples-1, :]
                #Trajectory information, traj::Array{Float32,2}, 1dim=DIM, 2dim=numsaples

                traj = ktraj[1:ndims, sample_range_start:sample_range_end]
                #Acquired data, data::Array{Complex{Float32},2}, 1dim=numsamples, 2dim=coils
                dat =  signal[sample_range_start:sample_range_end, :]

                #Saving profile
                push!(profiles, Profile(head, traj, dat))
            end
            #Update counters
            # scan_counter += 1
            scan_counter += adc_nSegments

            current += Nsamples
            if scan_counter % NadcsPerImage == 0 #For now only Nz is considered
                nz += 1 #another image
                scan_counter = 0 #reset counter
            end
        end
    end

    return RawAcquisitionData(params, profiles)
end