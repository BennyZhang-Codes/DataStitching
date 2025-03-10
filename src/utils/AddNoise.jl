
"""
    signal_noise = AddNoise(signal::AbstractArray{T, 2}, snr::Float64) where T <: Number
    Adds average white gaussian noise to the signal

# Arguments
- `signal::AbstractArray{T, 2} [nSample, nCha]`: signal data to be added with noise
- `snr::Float64`: signal to noise ratio

# Returns
- `signalout::AbstractArray{T, 2}`: signal with noise added

# Example
```julia
julia> signal = randn(520, 32)
julia> signal_noise = AddNoise(signal, 10.0)
```
"""
function AddNoise(signal::AbstractArray{T, 2}, snr::Real) where T <: Number
    nSample, nCha = size(signal);

    # Ensure input type consistency
    signalAmpl = sum(abs.(signal), dims=1) ./ nSample;

    # Target noise amplitude
    noiseAmpl = signalAmpl / snr

    if eltype(signal) <: Complex
        # Generate complex noise
        noise = (noiseAmpl / sqrt(2.0)) .* (randn(nSample, nCha) .+ 1im * randn(nSample, nCha))
    else
        # Generate real noise
        noise = noiseAmpl .* randn(nSample, nCha)
    end

    # Ensure the output matches the type of `signal`
    return signal .+ T.(noise)
end