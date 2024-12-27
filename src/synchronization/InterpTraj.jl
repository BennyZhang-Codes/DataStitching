"""
    traj_shifted = InterpTrajTime(traj::AbstractArray{<:Real, 2}, dt::Float64, delTime::Float64)
# Description
    Interpolates trajectory along the 1st dimension using Akima interpolation.

# Arguments
- `traj::AbstractArray{<:Real, 2}`: trajectory data (2D array, where each row is a point in the trajectory)
- `dt::Float64`: time delay between input samples
- `delTime::Float64`: time delay to be added or subtracted from time

# Keywords
- `intermode::Interpolations.InterpolationType`: interpolation mode (default is AkimaMonotonicInterpolation). 
        other Monotonic methods from Interpolations.jl are:
    1. LinearMonotonicInterpolation
    2. FiniteDifferenceMonotonicInterpolation
    3. CardinalMonotonicInterpolation
    4. AkimaMonotonicInterpolation
    5. FritschCarlsonMonotonicInterpolation
    6. FritschButlandMonotonicInterpolation
    7. SteffenMonotonicInterpolation

# Returns
- `traj_interp::Array{Float64, 2}`: interpolated trajectory data (2D array, where each row is a point in the interpolated trajectory)

# Examples
```julia-repl
julia> traj_shifted = InterpTrajTime(traj, dt, delTime)
```
"""
function InterpTrajTime(
    traj      ::AbstractArray{T, 2}      , 
    dt        ::Real                     , 
    delTime   ::Real                     ,
    datatime  ::AbstractVector{T}        ;
    intermode ::Interpolations.InterpolationType = AkimaMonotonicInterpolation()
    ) where {T<:Real}
    nSample, nTerm = size(traj)

    # Calculate adjusted time (delayed)
    trajTim = dt * (0:nSample-1) .- delTime

    # Perform Akima interpolation (assuming 1D interpolation for each row of the trajectory)
    traj_interp = []
    for i in 1:nTerm  # Loop over each column (dimension) of the trajectory
        itp = interpolate(trajTim, traj[:, i], intermode)
        etp = extrapolate(itp, 0)  # Extrapolate to the beginning of the time vector
        push!(traj_interp, etp.(datatime))  # Store the interpolated result
    end

    # Convert the result back to an array (if necessary)
    return T.(hcat(traj_interp...))  # Combine the interpolated columns into a single matrix
end

"""
    If dwell time is equal between trajectories and the sampling points.
    You can use this method.
"""
function InterpTrajTime(
    traj      ::AbstractArray{T, 2}      , 
    dt        ::Real                     , 
    delTime   ::Real                     ;
    intermode ::Interpolations.InterpolationType = AkimaMonotonicInterpolation()
    ) where {T<:Real}
    nSample, nTerm = size(traj)
    # Create time vector for input trajectory
    datatime = dt * (0:nSample-1)
    return InterpTrajTime(traj, dt, delTime, datatime; intermode=intermode)
end