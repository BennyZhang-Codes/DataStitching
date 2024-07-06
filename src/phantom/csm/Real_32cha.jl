


"""
```julia
>>> smap = csm_Real_32cha(200, 200)
>>> plot_imgs_subplots(abs.(smap), 4, 8)
```
"""
function csm_Real_32cha(Nx::Int64, Ny::Int64; verbose::Bool=false)
    if verbose
        @info "loading Real Coil-Sensitivity Map (32 channels)..."
    end
    sensitivity = MAT.matread("$(@__DIR__)/coilsensmap_32cha.mat")["coilsensmap_32cha"]
    out = imresize(sensitivity, (Nx, Ny, 32))
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
    out = out./ norm
    out[isnan.(out)] .= 0 + 0im
    return out
end