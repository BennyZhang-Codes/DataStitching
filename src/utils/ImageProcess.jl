# Normalization
function normalization(img::AbstractArray{T,2}) where T<:Real
    return (img.- minimum(img))./ (maximum(img) .- minimum(img))
end

# Standardization
function standardization(img::AbstractArray{T,2}) where T<:Real
    return (img.- mean(img))./ std(img; corrected=false)
end

"""
MSE
"""
function mymse(img::AbstractArray{T,2}, ref::AbstractArray{T,2}) where T<:Real
    return mean((img.- ref).^2)
end


function get_center_range(x::Int64, x_range::Int64)
    center = Int64(floor(x/2))
    return center - Int64(floor(x_range/2))+1 : center + Int64(ceil(x_range/2))
end


function  get_center_crop(images::Array{T, 3}, out_x::Int64, out_y::Int64) where T <: Number
    Nx = out_x
    Ny = out_y
    
    M, N, K = size(images)
    out = zeros(T, Nx, Ny, K)

    if Nx < M && Ny < N
        rangex = get_center_range(M, Nx)
        rangey = get_center_range(N, Ny)
        out = images[rangex, rangey, :]
    elseif Nx > M && Ny > N
        rangex = get_center_range(Nx, M)
        rangey = get_center_range(Ny, N)
        out[rangex, rangey, :] = images
    elseif Nx < M && Ny > N
        rangex = get_center_range(M, Nx)
        rangey = get_center_range(Ny, N)
        out[:, rangey, :] = images[rangex, :, :]
    elseif Nx > M && Ny < N
        rangex = get_center_range(Nx, M)
        rangey = get_center_range(N, Ny)
        out[rangex, :, :] = images[:, rangey, :]
    end
    return out
end
