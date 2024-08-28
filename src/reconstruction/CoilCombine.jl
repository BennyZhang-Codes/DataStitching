
function CoilCombineSOS(imgs::AbstractArray{T,N}, coildim::Int) where {T,N}
    # imgs: (nCoils, nZ, nY, nX)
    # coildim: 1
    imgs = sqrt.(sum(abs.(imgs).^2; dims=coildim))
    imgs = dropdims(imgs,dims=coildim)
    return imgs
end