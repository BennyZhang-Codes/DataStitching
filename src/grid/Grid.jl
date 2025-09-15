
export Grid
Base.@kwdef struct Grid{T<:AbstractFloat}
    nX::Int64 = 0
    nY::Int64 = 0
    nZ::Int64 = 0
    Δx::Real = 1.0
    Δy::Real = 1.0
    Δz::Real = 1.0
    x::AbstractVector{T} = [0.]
    y::AbstractVector{T} = [0.]
    z::AbstractVector{T} = [0.]
    matrixSize::Tuple{Int64,Int64,Int64} = (nX, nY, nZ)
    resolution::Tuple{Real,Real,Real} = (Δx, Δy, Δz)
end

@functor Grid

function Grid(
    nX::Int64, 
    nY::Int64, 
    nZ::Int64, 
    Δx::T, 
    Δy::T, 
    Δz::T;
    exchange_xy::Bool=false,
    reverse_x::Bool=false,
    reverse_y::Bool=false,
    reverse_z::Bool=false,
    ) where {T<:AbstractFloat}
    # x up->down, y left->right
    x, y, z = 1:nX, 1:nY, 1:nZ
    x, y = vec(x .+ y'*0.0), vec(x*0.0 .+ y') 
    x, y, z = vec(x .+ z'*0.0), vec(y .+ z'*0.0), vec(x*0.0 .+ z') #grid points
    x, y, z = x.-(nX+1)/2, y.-(nY+1)/2, z.-(nZ+1)/2
    x, y, z = x * Δx, y * Δy, z * Δz
    if exchange_xy
        x, y = y, x
    end
    if reverse_x
        x = reverse(x)
    end
    if reverse_y
        y = reverse(y)
    end
    if reverse_z
        z = reverse(z)
    end
    return Grid(nX=nX, nY=nY, nZ=nZ, Δx=Δx, Δy=Δy, Δz=Δz, x=T.(x), y=T.(y), z=T.(z))
end


Base.show(io::IO, b::Grid) = begin
    Δx = round(b.Δx*1e3, digits=2)
    Δy = round(b.Δy*1e3, digits=2)
    Δz = round(b.Δz*1e3, digits=2)
	print(io, "Grid [ MatrixSize: $(b.nX) x $(b.nY) x $(b.nZ), Resolution: $(Δx) x $(Δy) x $(Δz) mm³ ]")
end


export generate_grid
"""
    generate_grid(fov::Vector{<:Real}, matrix_size::Vector{<:Real}) -> Array{Float64,4}

Generate a 3D RPS (Read, Phase, Slice) grid of voxel center coordinates.

# Arguments
- `fov`: Field of view in millimeters, ordered as [FOV_R, FOV_P, FOV_S]
- `matrix_size`: Number of voxels along [R, P, S] axes

# Returns
- A 4D array of shape (nR, nP, nS, 3) where each voxel contains its RPS center coordinate
"""
function generate_grid(fov::Vector{<:Real}, matrix_size::Vector{<:Real})
    if length(fov) != 3 || length(matrix_size) != 3
        error("Both fov and matrix_size must be 3-element vectors.")
    end

    # Convert to Float64
    fov = Float64.(fov)
    matrix_size = Int.(matrix_size)

    spacing = fov ./ matrix_size

    r = ((0:matrix_size[1]-1) .- (matrix_size[1]-1)/2) .* spacing[1]
    p = ((0:matrix_size[2]-1) .- (matrix_size[2]-1)/2) .* spacing[2]
    s = ((0:matrix_size[3]-1) .- (matrix_size[3]-1)/2) .* spacing[3]

    rr = reshape(r, :, 1, 1)
    pp = reshape(p, 1, :, 1)
    ss = reshape(s, 1, 1, :)

    # Broadcast to form 3D grid
    grid = Array{Float64, 4}(undef, matrix_size[1], matrix_size[2], matrix_size[3], 3)
    grid[:, :, :, 1] .= rr
    grid[:, :, :, 2] .= pp
    grid[:, :, :, 3] .= ss

    return grid
end
