nX = 5
nY = 5
nZ = 3
Δx = Δy = Δz = 1.0

x, y, z = 1:nX, 1:nY, 1:nZ
# x up->down, y left->right
x, y, z = vec(x .+ y'*0.0), vec(x*0.0 .+ y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- nX/2 .- 1, y .- nY/2 .- 1



# x left->right, y down->up
x, y, z = vec(x*0.0 .+ y'), vec(x[end:-1:1] .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- (nX+1)/2, y .- (nY+1)/2


# x left->right, y up->down
x, y, z = vec(x*0.0 .+ y'), vec(x .+ 0.0*y'), vec(x*0.0 .+ y'*0.0) #grid points
x, y = x .- (nX+1)/2, y .- (nY+1)/2





x, y = vec(x .+ y'*0.0), vec(x*0.0 .+ y') 
x, y, z = vec(x .+ z'*0.0), vec(y .+ z'*0.0), vec(x*0.0 .+ z') #grid points
x, y, z = x.-(nX+1)/2, y.-(nY+1)/2, z.-(nZ+1)/2
x, y, z = x * Δx, y * Δy, z * Δz



Base.@kwdef struct Grid{T<:AbstractFloat}
    nX::Int64 = 1
    nY::Int64 = 1
    nZ::Int64 = 1
    Δx::T = 1.0
    Δy::T = 1.0
    Δz::T = 1.0
    x::AbstractVector{T} = [0.]
    y::AbstractVector{T} = [0.]
    z::AbstractVector{T} = [0.]
end

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


nX = 2
nY = 2
nZ = 1
Δx = Δy = Δz = 1.0

g = Grid(nX, nY, nZ, Δx, Δy, Δz)