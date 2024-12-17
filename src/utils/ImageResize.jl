"""
    imresize_real(img::AbstractArray{T}, dims::Tuple) where T<:Real
    Resize a real-valued image using imresize.
# Arguments
- `img::AbstractArray{T}`: The real-valued image to be resized.
- `dims::Tuple`: The new dimensions of the resized image.

# Examples
```julia
julia> img = rand(100, 100);
julia> plt_images(img, dim=2)
julia> plt_images(imresize_csm(img, (200, 200)), dim=2)
```
"""
function imresize_real(img::AbstractArray{T}, dims::Tuple) where T<:Real
    @assert length(dims) == ndims(img) "length of dims must be equal to ndims(img)"
    return imresize(img, dims...)
end

"""
    imresize_complex(img::AbstractArray{Complex{T}}, dims::Tuple) where T<:Real
    Resize a complex-valued image using imresize for the magnitude and angle separately.
# Arguments
- `img::AbstractArray{Complex{T}}`: The complex-valued image to be resized.
- `dims::Tuple`: The new dimensions of the resized image.

# Examples
```julia
julia> csm = randn(100, 100, 32) + 1im * randn(100, 100, 32);
julia> plt_images(abs.(csm), dim=3)
julia> plt_images(abs.(imresize_complex(csm, (200, 200, 32))), dim=3)
````
"""
function imresize_complex(img::AbstractArray{Complex{T}}, dims::Tuple) where T<:Real
    @assert length(dims) == ndims(img) "length of dims must be equal to ndims(img)"
    mag = imresize(abs.(img), dims...)
    pha = imresize(angle.(img), dims...)
    return mag .* exp.(im * pha)
end

