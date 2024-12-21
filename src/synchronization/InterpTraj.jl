"""
    traj_shifted = InterpTrajTime(traj::AbstractArray{<:Real, 2}, dt::Float64, delTime::Float64)
# Description
    Interpolates trajectory along the 1st dimension using Akima interpolation.

# Arguments
- `traj::AbstractArray{<:Real, 2}`: trajectory data (2D array, where each row is a point in the trajectory)
- `dt::Float64`: time delay between input samples
- `delTime::Float64`: time delay to be added or subtracted from time

# Returns
- `traj_interp::Array{Float64, 2}`: interpolated trajectory data (2D array, where each row is a point in the interpolated trajectory)

# Examples
```julia-repl
julia> traj_shifted = InterpTrajTime(traj, dt, delTime)
```
"""
function InterpTrajTime(
    traj::AbstractArray{<:Real, 2}, 
    dt::Float64, 
    delTime::Float64)
    nSample, nTerm = size(traj)

    # Create time vector for input trajectory
    datatime = dt * (0:nSample-1)

    # Calculate adjusted time (delayed)
    trajTim = datatime .- delTime

    # Perform Akima interpolation (assuming 1D interpolation for each row of the trajectory)
    traj_interp = []
    for i in 1:nTerm  # Loop over each column (dimension) of the trajectory
        itp = interpolate(trajTim, traj[:, i], AkimaMonotonicInterpolation())
        etp = extrapolate(itp, 0)  # Extrapolate to the beginning of the time vector
        push!(traj_interp, etp.(datatime))  # Store the interpolated result
    end

    # Convert the result back to an array (if necessary)
    return hcat(traj_interp...)  # Combine the interpolated columns into a single matrix
end