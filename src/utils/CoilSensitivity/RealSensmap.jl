


"""
```julia
>>> smap = RealCoilSensitivity_32cha(200, 200)
>>> plot_imgs_subplots(abs.(smap), 4, 8)
```
"""
function RealCoilSensitivity_32cha(Nx::Int64, Ny::Int64)
    sensitivity = MAT.matread("$(@__DIR__)/coilsensmap_32cha.mat")["coilsensmap_32cha"]
    out = imresize(sensitivity, (Nx, Ny, 32))
    norm = sqrt.(sum(abs.(out) .^ 2, dims=3))
    out = out./ norm
    return out
end