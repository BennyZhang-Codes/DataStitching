include("fan_mask.jl")
include("plot_imgs.jl")

export normalization, standardization

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